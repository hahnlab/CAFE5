#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <numeric>
#include <cmath>
#include <memory>

#include "utils.h" // for gene_family class
#include "clade.h"
#include "probability.h"
#include "matrix_cache.h"
#include "gene_family.h"

using namespace std;

extern std::mt19937 randomizer_engine;

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


#define GAMMA_CACHE_SIZE 1024

vector<double> lgamma_cache;

matrix chooseln_cache(100);

inline double lgamma2(double n)
{
    if (n >= 0 && n < GAMMA_CACHE_SIZE && (n - int(n) < 0.00000000001))
        return lgamma_cache.at(int(n));

    return lgamma(n);
}

void init_lgamma_cache()
{
    lgamma_cache.resize(GAMMA_CACHE_SIZE);
    for (int i = 0; i < GAMMA_CACHE_SIZE; ++i)
    {
        lgamma_cache[i] = lgamma(i);
    }

    for (int i = 0; i < chooseln_cache.size(); ++i)
        for (int j = 0; j < chooseln_cache.size(); ++j)
            chooseln_cache.set(i, j, lgamma2(i + 1) - lgamma2(j + 1) - lgamma2(i - j + 1));
}

double chooseln(double n, double r)
{
  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);

  if (n >= 0 && r >= 0 && n < chooseln_cache.size() && r < chooseln_cache.size() && (n - int(n) < 0.00000000001) && (r - int(r) < 0.00000000001))
      return chooseln_cache.get(n, r);

  return lgamma2(n + 1) - lgamma2(r + 1) - lgamma2(n - r + 1);
}

/* END: Math tools ----------------------- */

/* START: Birth-death model components --------------------- */

/* Eqn. (1) in 2005 paper. Assumes u = lambda 
  alpha should be lambda*t / 1+lambda*t
  coeff = 1 - 2 * alpha;
*/

#define MAX_FAMILY_SIZE_FOR_VECTORIZING 10000

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff)
{
    int m = std::min(c, s);
    double result = 0.0;
    if (m < MAX_FAMILY_SIZE_FOR_VECTORIZING)
    {
        int s_add_c = s + c;
        int s_add_c_sub_1 = s_add_c - 1;
        int s_sub_1 = s - 1;

        double t[MAX_FAMILY_SIZE_FOR_VECTORIZING];
        for (int j = 0; j <= m; j++) {
            t[j] = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
        }

        double expT[MAX_FAMILY_SIZE_FOR_VECTORIZING];
#ifdef HAVE_VECTOR_EXP
        vdExp(m + 1, t, expT);
#else
        std::transform(t, t + m + 1, expT, [](double d) { return exp(d); });
#endif

        double p[10000];
        for (int j = 0; j <= m; j++) {
            p[j] = expT[j] * pow(coeff, j); // Note that t is in log scale, therefore we need to do exp(t) to match Eqn. (1)
        }
        result = std::accumulate(p, p + m + 1, 0.0);
    }
    else
    {
        // Calculate sequentially rather than allocating the required memory to vectorize (much slower)
        double lastterm = 1;
        int s_add_c = s + c;
        int s_add_c_sub_1 = s_add_c - 1;
        int s_sub_1 = s - 1;

        for (int j = 0; j <= m; j++) {
            double t = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
            result += (exp(t) * lastterm); // Note that t is in log scale, therefore we need to do exp(t) to match Eqn. (1)
            lastterm *= coeff; // equivalent of ^j in Eqn. (1)
        }

    }
    return std::max(std::min(result, 1.0), 0.0);
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

//  if (result < .000000000000000001)
//  {
//    std::cout << "result= " << result << " lambda=" << lambda << " branch_length:" << branch_length << " From " << parent_size << " to " << size << std::endl;
//  }
  return result;
}


/* END: Birth-death model components ----------------------- */
void initialize_probabilities(const clade *p_tree, std::map<const clade *, std::vector<double> >& _probabilities, int _max_root_family_size, int _max_parsed_family_size)
{
    auto fn = [&](const clade *node) {
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

//! Calculates the probabilities of a given node for a given family size.
//! The probability of a leaf node is given from the family size at that node (1 for the family
//! size, and 0 for all other sizes). Internal nodes are calculated based on the probabilities
//! of the descendants.
void compute_node_probability(const clade *node, const gene_family&gene_family, const error_model*p_error_model,
    std::map<const clade *, std::vector<double> >& probabilities, 
    int _max_root_family_size,
    int _max_parsed_family_size,
    const lambda* _lambda,
    const matrix_cache& _calc) {
    if (node->is_leaf()) {
        int species_size = gene_family.get_species_size(node->get_taxon_name());

        if (p_error_model != NULL)
        {
            auto error_model_probabilities = p_error_model->get_probs(species_size);
            int offset = species_size - ((p_error_model->n_deviations() - 1) / 2);
            for (size_t i = 0; i<error_model_probabilities.size(); ++i)
            {
                if (offset + int(i) < 0)
                    continue;
                
                probabilities[node][offset+i] = error_model_probabilities[i];
            }
        }
        else
        {
            // cout << "Leaf node " << node->get_taxon_name() << " has " << _probabilities[node].size() << " probabilities" << endl;
            probabilities[node][species_size] = 1.0;
        }
    }

	else if (node->is_root()) {
		// at the root, the size of the vector holding the final likelihoods will be _max_root_family_size (size 0 is not included, so we do not add 1)
        std::vector<std::vector<double> > factors;
        auto fn = [&](const clade *c) {
            factors.push_back(_lambda->calculate_child_factor(_calc, c, probabilities[c], 1, _max_root_family_size, 0, _max_parsed_family_size));
        };
        node->apply_to_descendants(fn);
        vector<double>& node_probs = probabilities[node];
        // factors[0] is left child
        // factors[1] is right child
        for (size_t i = 0; i < node_probs.size(); i++) {
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
        auto fn = [&](const clade *c) {
            factors.push_back(_lambda->calculate_child_factor(_calc, c, probabilities[c], 0, _max_parsed_family_size, 0, _max_parsed_family_size));
        };

        node->apply_to_descendants(fn);

        vector<double>& node_probs = probabilities[node];
        // factors[0] is left child
        // factors[1] is right child
        for (size_t i = 0; i < node_probs.size(); i++) {
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

std::vector<double> get_random_probabilities(const clade *p_tree, int number_of_simulations, int root_family_size, int max_family_size, int max_root_family_size, const lambda *p_lambda, const matrix_cache& cache, error_model *p_error_model)
{
    vector<double> result(number_of_simulations);
    vector<gene_family> families(number_of_simulations);
    vector<std::map<const clade *, std::vector<double> >> pruners(number_of_simulations);

    // TODO: This is slow so it should be done in parallel. Care will have to be taken
    // that stl containers that are added to are thread-safe.
    for (size_t i = 0; i < result.size(); ++i)
    {
        clademap<int> sizes;
        auto fn = [&sizes, p_lambda, p_error_model, max_family_size, &cache](const clade *c)
        {
            set_weighted_random_family_size(c, &sizes, p_lambda, p_error_model, max_family_size, cache);
        };

        sizes[p_tree] = root_family_size;
        p_tree->apply_prefix_order(fn); 

        for (auto& it : sizes) {
            if (it.first->is_leaf())
            {
                families[i].set_species_size(it.first->get_taxon_name(), it.second);
            }
        }
        // pruners[i].reset(new likelihood_computer(max_root_family_size, max_family_size, p_lambda, families[i], cache, NULL));
        initialize_probabilities(p_tree, pruners[i], max_root_family_size, max_family_size);
    }

#pragma omp parallel for
    for (size_t i = 0; i < result.size(); ++i)
    {
        auto fn = [&](const clade *c) { compute_node_probability(c, families[i], NULL, pruners[i], max_root_family_size, max_family_size, p_lambda, cache); };
        p_tree->apply_reverse_level_order(fn);
        result[i] = *std::max_element(pruners[i].at(p_tree).begin(), pruners[i].at(p_tree).end());
    }

    sort(result.begin(), result.end());

    return result;
}

//! Set the family size of a node to a random value, using parent's family size
void set_weighted_random_family_size(const clade *node, clademap<int> *sizemap, const lambda *p_lambda, error_model *p_error_model, int max_family_size, const matrix_cache& cache)
{
    if (node->is_root()) // if node is root, we do nothing
        return;

    int parent_family_size = (*sizemap)[node->get_parent()];
    size_t c = 0; // c is the family size we will go to

    double lambda = p_lambda->get_value_for_clade(node);
    double branch_length = node->get_branch_length();

    if (parent_family_size > 0) {
        auto probabilities = cache.get_matrix(branch_length, lambda);
        if (cache.is_saturated(branch_length, lambda))
        {
            std::uniform_int_distribution<int> distribution(0, max_family_size - 1);
            c = distribution(randomizer_engine);
        }
        vector<double> v;
        for (int i = 0; i < max_family_size; i++) {
            v.push_back(probabilities->get(parent_family_size, i));
        }
        std::discrete_distribution<int> distribution(v.begin(), v.end());
        c = distribution(randomizer_engine);
    }

    if (node->is_leaf())
    {
        c = adjust_for_error_model(c, p_error_model);
    }

    (*sizemap)[node] = c;
}

size_t adjust_for_error_model(size_t c, const error_model *p_error_model)
{
    if (p_error_model == nullptr)
        return c;

    if (c >= p_error_model->get_max_count())
    {
        throw runtime_error("Trying to simulate leaf family size that was not included in error model");
    }
    auto probs = p_error_model->get_probs(c);

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rnd = distribution(randomizer_engine);
    if (rnd < probs[0])
    {
        c--;
    }
    else if (rnd >(1 - probs[2]))
    {
        c++;
    }

    return c;
}
