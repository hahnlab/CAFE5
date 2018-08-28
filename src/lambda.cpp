#include <vector>
#include <map>
#include <cmath>
#include <iomanip>

#include "lambda.h"
#include "clade.h"
#include "fminsearch.h"
#include "utils.h"
#include "probability.h"
#include "core.h"
#include "matrix_cache.h"

using namespace std;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

std::vector<double> single_lambda::calculate_child_factor(const matrix_cache& calc, clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const
{
    auto matrix = calc.get_matrix(probabilities.size(), child->get_branch_length(), _lambda);
#if 0
    printf("  Node %s matrix parameters: %d, %f, %f\n", child->get_taxon_name().c_str(), probabilities.size(), child->get_branch_length(), _lambda);
    printf("  Multipliers: %d, %d, %d, %d\n", s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);
    if (matrix.size() > 65)
        printf("  Matrix 65, 65 is: %e\n", matrix.get(65, 65));
#endif
	return matrix.multiply(probabilities, s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);
}

std::string single_lambda::to_string() const
{
    ostringstream ost;
    ost << setw(15) << setprecision(14) << _lambda;
    return ost.str();
}

std::vector<double> multiple_lambda::calculate_child_factor(const matrix_cache& calc, clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const
{
	std::string nodename = child->get_taxon_name();
	int lambda_index = _node_name_to_lambda_index.at(nodename);
	double lambda = _lambdas[lambda_index];
	//cout << "Matrix for " << child->get_taxon_name() << endl;
	auto matrix = calc.get_matrix(probabilities.size(), child->get_branch_length(), lambda); // Ben: is _factors[child].size() the same as _max_root_family_size? If so, why not use _max_root_family_size instead?
	return matrix.multiply(probabilities, s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);
}

void multiple_lambda::update(double* values)
{
    std::copy(values, values + _lambdas.size(), _lambdas.begin());
}

std::string multiple_lambda::to_string() const
{
    ostringstream ost;
    ost << setw(15) << setprecision(14);
    for (size_t i = 0; i < _lambdas.size(); ++i)
    {
        ost << _lambdas[i];
        if (i != _lambdas.size() - 1) ost << ", ";
    }
    return ost.str();
}

bool multiple_lambda::is_valid()
{
    return std::none_of(_lambdas.begin(), _lambdas.end(), [](double d) { return d < 0; });
}

double multiple_lambda::get_value_for_clade(clade *c) const {
    int index = _node_name_to_lambda_index.at(c->get_taxon_name());
    return _lambdas[index];
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

double fn_calc_score(double* p_lambda, void* args)
{
    optimizer_scorer *opt = reinterpret_cast<optimizer_scorer*>(args);
    opt->calculate_score(p_lambda);
}

optimizer::optimizer(optimizer_scorer *p_scorer) : _p_scorer(p_scorer)
{
#ifdef SILENT
    quiet = true;
#endif
    pfm = fminsearch_new();
}

optimizer::result optimizer::optimize()
{
    auto initial = _p_scorer->initial_guesses();
    fminsearch_set_equation(pfm, fn_calc_score, initial.size(), _p_scorer);
    if (explode)
    {
        pfm->rho = 1.5;				// reflection
        pfm->chi = 50;				// expansion
        pfm->delta = 0.4;
    }
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
    pfm->maxiters = 25;
	fminsearch_min(pfm, &initial[0]);
    double *re = fminsearch_get_minX(pfm);

    result r;
    r.score = fminsearch_get_minF(pfm);
    r.values.resize(initial.size());
    copy(re, re + initial.size(), r.values.begin());
    if (!quiet)
    {
        log_results(pfm, initial, re);
    }
    return r;
}

void optimizer::log_results(FMinSearch * pfm, std::vector<double> &initial, double * re)
{
    if (fminsearch_get_minF(pfm) == -log(0))
    {
        cerr << "Failed to find any reasonable values" << endl;
    }
    else
    {
        cout << "Completed " << pfm->iters << " iterations" << endl;
        cout << "Best match" << (initial.size() == 1 ? " is: " : "es are: ") << setw(15) << setprecision(14);
        for (size_t i = 0; i < initial.size(); ++i)
            cout << re[i] << ',';
        cout << endl;
    }
}


