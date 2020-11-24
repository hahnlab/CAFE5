#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <numeric>
#include <cmath>
#include <memory>

#include "easylogging++.h"


#ifdef HAVE_VECTOR_EXP
#include "mkl.h"
#endif

#include "clade.h"
#include "probability.h"
#include "matrix_cache.h"
#include "gene_family.h"
#include "error_model.h"

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

/// <summary>
/// Call this method if there are less than MAX_POSSIBLE_FAMILY_SIZE sizes. Uses strictly stack memory 
/// to run faster and in a thread-compatible manner
void compute_node_probability_small_families(const clade *node, const gene_family&gene_family, const error_model*p_error_model,
    std::map<const clade *, std::vector<double> >& probabilities, 
    int max_root_family_size,
    int max_family_size,
    const lambda* _lambda,
    const matrix_cache& _calc) {

    assert(max(max_root_family_size, max_family_size) < MAX_STACK_FAMILY_SIZE);

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
        vector<double>& node_probs = probabilities[node];
        fill(node_probs.begin(), node_probs.end(), 1);

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            double result[MAX_STACK_FAMILY_SIZE];
             _lambda->calculate_child_factor(_calc, *it, probabilities[*it], 1, max_root_family_size, 0, max_family_size, result);
             for (size_t i = 0; i < node_probs.size(); i++) {
                 node_probs[i] *= result[i];
             }
        }
    }
    else {
        // at any internal node, the size of the vector holding likelihoods will be _max_parsed_family_size+1 because size=0 is included
        vector<double>& node_probs = probabilities[node];
        fill(node_probs.begin(), node_probs.end(), 1);

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            double result[MAX_STACK_FAMILY_SIZE];
            _lambda->calculate_child_factor(_calc, *it, probabilities[*it], 0, max_family_size, 0, max_family_size, result);
            for (size_t i = 0; i< node_probs.size(); i++) {
                node_probs[i] *= result[i];
            }
        }
    }
}

/// <summary>
/// Call this method if there are more than MAX_POSSIBLE_FAMILY_SIZE sizes. Allocates sufficient heap memory 
/// to handle ay number of families
void compute_node_probability_large_families(const clade* node, const gene_family& gene_family, const error_model* p_error_model,
    std::map<const clade*, std::vector<double> >& probabilities,
    int max_root_family_size,
    int max_family_size,
    const lambda* _lambda,
    const matrix_cache& _calc) {
    if (node->is_leaf()) {
        int species_size = gene_family.get_species_size(node->get_taxon_name());

        if (p_error_model != NULL)
        {
            auto error_model_probabilities = p_error_model->get_probs(species_size);
            int offset = species_size - ((p_error_model->n_deviations() - 1) / 2);
            for (size_t i = 0; i < error_model_probabilities.size(); ++i)
            {
                if (offset + int(i) < 0)
                    continue;

                probabilities[node][offset + i] = error_model_probabilities[i];
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
        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            vector<double> result(max_root_family_size);
            _lambda->calculate_child_factor(_calc, *it, probabilities[*it], 1, max_root_family_size, 0, max_family_size, result.data());
            factors.push_back(result);
        }
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

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            vector<double> result(max_family_size+1);
            _lambda->calculate_child_factor(_calc, *it, probabilities[*it], 0, max_family_size, 0, max_family_size, result.data());
            factors.push_back(result);
        }

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

//! Calculates the probabilities of a given node for a given family size.
//! The probability of a leaf node is given from the family size at that node (1 for the family
//! size, and 0 for all other sizes). Internal nodes are calculated based on the probabilities
//! of the descendants.
void compute_node_probability(const clade* node, const gene_family& gene_family, const error_model* p_error_model,
    std::map<const clade*, std::vector<double> >& probabilities,
    int max_root_family_size,
    int max_family_size,
    const lambda* _lambda,
    const matrix_cache& calc) {
    if (max(max_root_family_size, max_family_size) < MAX_STACK_FAMILY_SIZE)
        compute_node_probability_small_families(node, gene_family, p_error_model, probabilities, max_root_family_size, max_family_size, _lambda, calc);
    else
        compute_node_probability_large_families(node, gene_family, p_error_model, probabilities, max_root_family_size, max_family_size, _lambda, calc);
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

/// <summary>
/// Create a gene family based on the given tree and root size, by assigning child sizes down the tree 
/// The child size is selected randomly, based on the lambda and error model
/// </summary>
gene_family create_family(const clade* p_tree, int root_family_size, int max_family_size, const lambda* p_lambda, const matrix_cache& cache)
{
    gene_family result;
    // generate a tree with root_family_size at the root
    clademap<int> sizes;
    sizes[p_tree] = root_family_size;

    // note we do not use an error model for creating family sizes. See architecture decision #6
    p_tree->apply_prefix_order([&](const clade* c) { set_weighted_random_family_size(c, &sizes, p_lambda, nullptr, max_family_size, cache); });
    result.init_from_clademap(sizes);
    return result;
}

vector<double> compute_family_probabilities(const clade* p_tree, const vector<gene_family>& families, int max_family_size, int max_root_family_size, const lambda* p_lambda, const matrix_cache& cache)
{
    vector<double> result(families.size());

    // Allocate space to calculate all of the families simultaneously
    vector<clademap<std::vector<double>>> pruners(families.size());
    for (auto& p : pruners)
    {
        // vector of lk's at tips must go from 0 -> _max_possible_family_size, so we must add 1
        auto fn = [&](const clade* node) { p[node].resize(node->is_root() ? max_root_family_size : max_family_size + 1); };
        for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), fn);
    }

#pragma omp parallel for
    for (size_t i = 0; i < result.size(); ++i)
    {
        for (auto it = p_tree->reverse_level_begin(); it != p_tree->reverse_level_end(); ++it)
            compute_node_probability(*it, families[i], NULL, pruners[i], max_root_family_size, max_family_size, p_lambda, cache);
        result[i] = *std::max_element(pruners[i].at(p_tree).begin(), pruners[i].at(p_tree).end());
    }
    return result;
}

/*! Create a sorted vector of probabilities by generating random trees 
    \param p_tree The structure of the tree to generate
    \param number_of_simulations The number of random probabilities to return
    \param root_family_size The count of the family at the root. All other family sizes will be generated randomly based on this
    \param max_family_size The maximum possible family size (Used to cut off probability calculations at a reasonable value)
    \param root_family_size The maximum possible family size at the root
    \param lambda The rate of change of the family

    \returns a sorted vector of probabilities of the requested number of randomly generated trees
*/

std::vector<double> get_random_probabilities(const clade *p_tree, int number_of_simulations, int root_family_size, int max_family_size, int max_root_family_size, const lambda *p_lambda, const matrix_cache& cache)
{
    VLOG(2) << "Getting random probabilities for " << number_of_simulations << " simulations";
    vector<gene_family> families(number_of_simulations);

    generate(families.begin(), families.end(), [&]() { return create_family(p_tree, root_family_size, max_family_size, p_lambda, cache); });

    auto result = compute_family_probabilities(p_tree, families, max_family_size, max_root_family_size, p_lambda, cache);

    sort(result.begin(), result.end());

    VLOG(2) << "Random probabilities complete";

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
        c = probabilities->select_random_y(parent_family_size, max_family_size);
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

    if (c >= p_error_model->get_max_family_size())
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

double pvalue(double v, const vector<double>& conddist)
{
    int idx = conddist.size() - 1;

    auto bound = std::upper_bound(conddist.begin(), conddist.end(), v);
    if (bound != conddist.end())
    {
        idx = bound - conddist.begin();
    }
    return  idx / (double)conddist.size();
}

//! Compute pvalues for each family based on the given lambda
vector<double> compute_pvalues(const clade* p_tree, const std::vector<gene_family>& families, const lambda* p_lambda, const matrix_cache& cache, int number_of_simulations, int max_family_size, int max_root_family_size)
{
    LOG(INFO) << "Computing pvalues...";

    const int mx = max_family_size;
    const int mxr = max_root_family_size;

    std::vector<std::vector<double> > conditional_distribution(mxr);
    for (int i = 0; i < mxr; ++i)
    {
        conditional_distribution[i] = get_random_probabilities(p_tree, number_of_simulations, i, mx, mxr, p_lambda, cache);
    }
    VLOG(1) << "Conditional distributions calculated";

    auto observed_max_likelihoods = compute_family_probabilities(p_tree, families, max_family_size, max_root_family_size, p_lambda, cache);

    vector<double> result(observed_max_likelihoods.size());

    transform(observed_max_likelihoods.begin(), observed_max_likelihoods.end(), result.begin(), [&](double observed_max_likelihood)
        {
            vector<double> pvalues(mxr);
            for (size_t s = 0; s < (size_t)mxr; s++)
            {
                pvalues[s] = pvalue(observed_max_likelihood, conditional_distribution[s]);
            }

            return *max_element(pvalues.begin() + 1, pvalues.end());
        });


    LOG(INFO) << "done!\n";

    return result;
}
