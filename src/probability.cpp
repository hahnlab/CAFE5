#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>

#include "utils.h" // for gene_family class
#include "clade.h"
#include "probability.h"
#include "matrix_cache.h"
#include "gene_family.h"

/* Useful links
1) http://www.rskey.org/gamma.htm # explanation for lgamma
2) http://www.physics.unlv.edu/~pang/cp_c.html # c code
*/

/* Necessary for old C implementation of what now is lgamma
#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
*/

/* Necessary for old C implementation of what now is lgamma 
double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
 24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
 -5.395239384953e-6 };
*/

/* Old C implementation of what now is lgamma */
/*
double gammaln(double a)
{
	int n;
	double p = __Qs[0];
	double a_add_5p5 = a + 5.5;
	for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
	return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}
*/


/* START: Math tools --------------------- */
/* Old C implementation necessary for set_node_familysize_random. Now using uniform_real_distribution()
*/
double unifrnd() {
  double result = rand() / (RAND_MAX + 1.0); // rand() returns an int from 0 to RAND_MAX (which is defined in std); the +1.0 is there probably so that we do not draw exactly 1.
  return result;
}

vector<double> lgamma_cache;
void init_lgamma_cache()
{
    lgamma_cache.resize(512);
    for (int i = 0; i < 512; ++i)
    {
        lgamma_cache[i] = lgamma(i);
    }
}

inline double lgamma2(double n)
{
    if (n < 512 && (n - int(n) < 0.00000000001))
        return lgamma_cache.at(int(n));

    return lgamma(n);
}

double chooseln(double n, double r)
{
  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);
  return lgamma2(n + 1) - lgamma2(r + 1) - lgamma2(n - r + 1);
}

/* END: Math tools ----------------------- */

/* START: Birth-death model components --------------------- */

/* Eqn. (1) in 2005 paper. Assumes u = lambda 
  alpha should be lambda*t / 1+lambda*t
  coeff = 1 - 2 * alpha;
*/

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff)
{
    int m = std::min(c, s);
    double lastterm = 1;
    double p = 0.0;
    int s_add_c = s + c;
    int s_add_c_sub_1 = s_add_c - 1;
    int s_sub_1 = s - 1;

    for (int j = 0; j <= m; j++) {
        double t = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
        p += (exp(t) * lastterm); // Note that t is in log scale, therefore we need to do exp(t) to match Eqn. (1)
        lastterm *= coeff; // equivalent of ^j in Eqn. (1)
    }

    return std::max(std::min(p, 1.0), 0.0);
  return std::max(std::min(p, 1.0), 0.0);
}

double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size)
{
  //std::cout << "Lambda: " << lambda << ", branch_length: " << branch_length;
  double alpha = lambda*branch_length / (1 + lambda*branch_length);
  double coeff = 1 - 2 * alpha;

  double result = 0;
  if (coeff > 0 && coeff != 1)  // if coeff < 0 branch is saturated (some of the characters may have changed and then changed again)
  {
      result = birthdeath_rate_with_log_alpha(parent_size, size, log(alpha), coeff);
  }
      
  //printf("Birthdeath rate for 1, 0 (alpha=%f, coeff=%f), : %f\n", alpha, coeff, birthdeath_rate_with_log_alpha(1, 0, log(alpha), coeff));


//  if (result < .000000000000000001)
//  {
//    std::cout << "result= " << result << " lambda=" << lambda << " branch_length:" << branch_length << " From " << parent_size << " to " << size << std::endl;
//  }
  return result;
}

/* END: Birth-death model components ----------------------- */

void likelihood_computer::initialize_memory(const clade *p_tree)
{
    auto fn = [this](const clade *node) {
        if (node->is_root())
        {
            _probabilities[node].resize(_max_root_family_size);
        }
        else
        {
            _probabilities[node].resize(_max_parsed_family_size + 1); // vector of lk's at tips must go from 0 -> _max_possible_family_size, so we must add 1
        }
    };
    p_tree->apply_reverse_level_order(fn);
}

//! Operator () overload for likelihood_computer.
/*!
  The operator () overload here allows likelihood_computer to be called as a function (i.e., it makes likelihood_computer a functor).
  This is what allows likelihood_computer to be called recursively in the pruning algorithm, through apply_reverse_level_order.
*/
void likelihood_computer::operator()(const clade *node) {
    if (node->is_leaf()) {
        int species_size = _gene_family.get_species_size(node->get_taxon_name());

        if (_p_error_model != NULL)
        {
            auto error_model_probabilities = _p_error_model->get_probs(species_size);
            int offset = species_size - ((_p_error_model->n_deviations() - 1) / 2);
            for (size_t i = 0; i<error_model_probabilities.size(); ++i)
            {
                if (offset + int(i) < 0)
                    continue;
                
                _probabilities[node][offset+i] = error_model_probabilities[i];
            }
        }
        else
        {
            // cout << "Leaf node " << node->get_taxon_name() << " has " << _probabilities[node].size() << " probabilities" << endl;
            _probabilities[node][species_size] = 1.0;
        }
    }

	else if (node->is_root()) {
		// at the root, the size of the vector holding the final likelihoods will be _max_root_family_size (size 0 is not included, so we do not add 1)
        std::vector<std::vector<double> > factors;
        auto fn = [this, &factors](const clade *c) {
            factors.push_back(_lambda->calculate_child_factor(_calc, c, _probabilities[c], 1, _max_root_family_size, 0, _max_parsed_family_size));
        };
        node->apply_to_descendants(fn);
        vector<double>& node_probs = _probabilities[node];
        // factors[0] is left child
        // factors[1] is right child
        double s = 1.0;
        for (int i = 0; i < node_probs.size(); i++) {
            node_probs[i] = 1;
            auto it = factors.begin();

            for (; it != factors.end(); it++) {
                node_probs[i] *= it->at(i);
            }
        }
    }

    else {
		// at any internal node, the size of the vector holding likelihoods will be _max_parsed_family_size+1 because size=0 is included
        std::vector<std::vector<double> > factors;
        auto fn = [this, &factors](const clade *c) {
            factors.push_back(_lambda->calculate_child_factor(_calc, c, _probabilities[c], 0, _max_parsed_family_size, 0, _max_parsed_family_size));
        };

        node->apply_to_descendants(fn);

        vector<double>& node_probs = _probabilities[node];
        // factors[0] is left child
        // factors[1] is right child
        double s = 1.0;
        for (int i = 0; i < node_probs.size(); i++) {
            node_probs[i] = 1;
            auto it = factors.begin();

            for (; it != factors.end(); it++) {
                node_probs[i] *= it->at(i);
            }
        }
    }
}

/* END: Likelihood computation ---------------------- */

std::vector<int> uniform_dist(int n_draws, int min, int max) {
    
    std::random_device rd; // seed
    std::default_random_engine generator(rd()); // seeding generator
    std::uniform_int_distribution<int> distribution(min, max); // initializing uniform generator
    std::vector<int> uniform_vec(n_draws); // for storing results
    
    for (int i = 0; i < n_draws; ++i) {
        int number = distribution(generator); // drawing from uniform generator by plugging in random number
        //cout << "Number is: " << number << endl;
        uniform_vec[i] = number;
    }
        
    return uniform_vec;
}

/* START: Weighted draw from vector */
//! Draw ints or doubles from n_draws equal intervals using specified weights.
void category_selector::weighted_cat_draw(int n_draws, std::vector<double> gamma_cat_probs) {
    std::random_device rd;
    std::mt19937 gen(rd()); // seeding random number engine

    // creating equal-sized intervals (n_draws of them)
    std::vector<double> intervals(gamma_cat_probs.size() + 1);
    for (int i = 0; i != intervals.size(); ++i) {
        intervals[i] = i;
        //std::cout << i << " ith interval" << std::endl;
    }

    std::piecewise_constant_distribution<double> d(intervals.begin(), intervals.end(), gamma_cat_probs.begin());

    _gamma_cats.resize(n_draws);
    // now drawing
    for (int i = 0; i < n_draws; ++i) {
        _gamma_cats[i] = d(gen);
        //cout << (*p_gamma_cats)[i] << ", ";
    }
}
/* END: Weighted draw from vector */

std::vector<double> get_random_probabilities(const clade *p_tree, int number_of_simulations, int root_family_size, int max_family_size, const lambda *p_lambda, const matrix_cache& cache, error_model *p_error_model)
{
    vector<double> result(number_of_simulations);
    vector<gene_family> families(number_of_simulations);
    vector<unique_ptr<likelihood_computer>> pruners(number_of_simulations);

    // TODO: This is slow so it should be done in parallel. Care will have to be taken
    // that stl containers that are added to are thread-safe.
    for (size_t i = 0; i < result.size(); ++i)
    {
        clademap<int> sizes;
        random_familysize_setter rfs(&sizes, max_family_size, p_lambda, p_error_model, cache);
        sizes[p_tree] = root_family_size;
        p_tree->apply_prefix_order(rfs); // this is where the () overload of random_familysize_setter is used

        for (auto& it : sizes) {
            if (it.first->is_leaf())
            {
                families[i].set_species_size(it.first->get_taxon_name(), it.second);
            }
        }
        pruners[i].reset(new likelihood_computer(root_family_size, max_family_size, p_lambda, families[i], cache, NULL));
        pruners[i]->initialize_memory(p_tree);
    }

#pragma omp parallel for
    for (size_t i = 0; i < result.size(); ++i)
    {
        p_tree->apply_reverse_level_order(*pruners[i]);
        result[i] = pruners[i]->max_likelihood(p_tree);
    }

    sort(result.begin(), result.end());

    return result;
}

std::vector<std::vector<double> > get_conditional_distribution_matrix(const clade *p_tree, int root_family_size, int max_family_size, int number_of_simulations, const lambda * p_lambda, const matrix_cache& cache)
{
    // sanity check to verify the matrix multiplications won't cause memory corruption
    int s = cache.get_matrix_size();
    if (s < root_family_size || s < max_family_size)
        throw runtime_error("Precalculated matrices too small");
    
    std::vector<std::vector<double> > matrix(root_family_size);
	for (int i = 0; i < root_family_size; ++i)
	{
		matrix[i] = get_random_probabilities(p_tree, number_of_simulations, root_family_size, max_family_size, p_lambda, cache, NULL);
	}
	return matrix;
}
