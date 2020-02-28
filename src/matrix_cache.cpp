#include <algorithm>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <random>

#include "matrix_cache.h"
#include "probability.h"

#ifdef HAVE_BLAS
#include "mkl.h"
#endif

extern std::mt19937 randomizer_engine;

bool matrix::is_zero() const
{
    return *max_element(values.begin(), values.end()) == 0;
}

//! Take in a matrix and a vector, compute product, return it
/*!
This function returns a likelihood vector by multiplying an initial likelihood vector and a transition probability matrix.
A minimum and maximum on the parent's and child's family sizes is provided. Because the root is forced to be >=1, for example, s_min_family_size for the root could be set to 1.
*/
vector<double> matrix::multiply(const vector<double>& v, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const
{
    vector<double> result(s_max_family_size - s_min_family_size + 1);

    assert(c_min_family_size < c_max_family_size);
    assert(v.size() > size_t(c_max_family_size - c_min_family_size));

#ifdef HAVE_BLAS
    double alpha = 1.0, beta = 0.;
    int m = s_max_family_size - s_min_family_size + 1;
    int k = c_max_family_size - c_min_family_size + 1;
    int n = 1;
    const double *sub = &values[0] + s_min_family_size*_size + c_min_family_size;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, sub, _size, &v[0], n, beta, &result[0], n);
#else
    //cout << "Matrix multiply " << matrix.size() << "x" << v.size() << " (submatrix " << s_min_family_size << ":" << s_max_family_size;
    //cout << " " << c_min_family_size << ":" << c_max_family_size << ")" << endl;

    //assert(s_max_family_size - s_min_family_size == c_max_family_size - c_min_family_size);

    for (int s = s_min_family_size; s <= s_max_family_size; s++) {
        result[s - s_min_family_size] = 0;

        for (int c = c_min_family_size; c <= c_max_family_size; c++) {
            result[s - s_min_family_size] += get(s, c) * v[c - c_min_family_size];
        }
    }
#endif
    return result;
}

int matrix::select_random_y(int x, int max) const
{
    assert(x < _size);
    assert(max < _size);
    std::discrete_distribution<int> distribution(values.begin() + x*_size, values.begin() + x*_size + max);
    return distribution(randomizer_engine);
}

matrix_cache::~matrix_cache()
{
    for (auto m : _matrix_cache)
    {
        delete m.second;
    }
}

/* START: Likelihood computation ---------------------- */
//! Calls BD formulas, but checks/populates cache
double matrix_cache::get_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int child_size) const {
    // The probability of 0 remaining 0 is 1
    // The probability of 0 going to any other count is 0 (if you lose the gene family, you do not regain it)
    if (parent_size == 0.0)
        return child_size == 0.0 ? 1.0 : 0.0;

    return the_probability_of_going_from_parent_fam_size_to_c(lambda, branch_length, parent_size, child_size);
}

//! Compute transition probability matrix for all gene family sizes from 0 to size-1 (=_max_root_family_size-1)
const matrix* matrix_cache::get_matrix(double branch_length, double lambda) const {
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    matrix *result = NULL;
    matrix_cache_key key(_matrix_size, lambda, branch_length);
    if (_matrix_cache.find(key) != _matrix_cache.end())
    {
        result = _matrix_cache.at(key);
    }

    if (result == NULL)
    {
        ostringstream ost;
        ost << "Failed to find matrix for " << _matrix_size << "," << branch_length << "," << lambda;
        throw std::runtime_error(ost.str());
    }
    return result;
}

vector<double> get_lambda_values(const lambda *p_lambda)
{
    vector<double> lambdas;
    auto sl = dynamic_cast<const single_lambda *>(p_lambda);
    if (sl)
    {
        lambdas.push_back(sl->get_single_lambda());
    }
    else
    {
        auto ml = dynamic_cast<const multiple_lambda *>(p_lambda);
        lambdas = ml->get_lambdas();
    }
    return lambdas;
}

bool matrix_cache::is_saturated(double branch_length, double lambda)
{
    double alpha = lambda*branch_length / (1 + lambda*branch_length);
    return (1 - 2 * alpha) < 0;
}

void matrix_cache::precalculate_matrices(const std::vector<double>& lambdas, const std::set<double>& branch_lengths)
{
	// build a list of required matrices
	vector<matrix_cache_key> keys;
	for (double lambda : lambdas)
	{
		for (double branch_length : branch_lengths)
		{
			matrix_cache_key key(_matrix_size, lambda, branch_length);
			if (_matrix_cache.find(key) == _matrix_cache.end())
			{
				keys.push_back(key);
			}
		}
	}

	// calculate matrices in parallel
	vector<matrix*> matrices(keys.size());
	generate(matrices.begin(), matrices.end(), [this] { return new matrix(this->_matrix_size); });

	int s = 0;
	size_t i = 0;
	size_t num_keys = keys.size();

    auto t1 = std::chrono::high_resolution_clock::now();

#pragma omp parallel for private(s) collapse(2)
	for (i = 0; i < num_keys; ++i)
	{
		for (s = 1; s < _matrix_size; s++) {
			double lambda = keys[i].lambda();
			double branch_length = keys[i].branch_length();

			matrix* m = matrices[i];
			m->set(0, 0, get_from_parent_fam_size_to_c(lambda, branch_length, 0, 0));
			if (!is_saturated(branch_length, lambda))
			{
				for (int j = 0; j < m->size(); ++j)
				{
					m->set(0, j, get_from_parent_fam_size_to_c(lambda, branch_length, 0, j));
				}
				for (int c = 0; c < _matrix_size; c++) {
					m->set(s, c, get_from_parent_fam_size_to_c(lambda, branch_length, s, c));
				}
			}
		}
	}

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "matrix precalcuation time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;

    // copy matrices to our internal map
    for (size_t i = 0; i < keys.size(); ++i)
    {
        _matrix_cache[keys[i]] = matrices[i];
    }
}

void matrix_cache::warn_on_saturation(std::ostream& ost)
{
    for (auto& kv : _matrix_cache)
    {
        if (is_saturated(kv.first.branch_length(), kv.first.lambda()))
            ost << "WARNING: Saturated branch using lambda " << kv.first.lambda() << " on branch length " << kv.first.branch_length() << endl;
    }
}

std::ostream& operator<<(std::ostream& ost, matrix_cache& c)
{
    ost << c.get_cache_size() << " matrices. Keys: ";
    for (auto& kv : c._matrix_cache)
    {
        ost << "(" << kv.first.branch_length() << "," << kv.first.lambda() << "),";
    }
    return ost;
}

