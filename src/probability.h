#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <tuple>
#include <set>

#include <iostream>
#include <set>
#include "io.h"
#include "lambda.h"
#include "clade.h"

class matrix;
class matrix_cache;

struct pvalue_parameters
{
    const clade* p_tree;
    const lambda* p_lambda;
    const int max_family_size;
    const int max_root_family_size;
    const matrix_cache& cache;
};


double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff);
double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size);
double chooseln(double n, double k);

/* START: Likelihood computation ---------------------- */

/// compute_node_probability is a function that is used in the pruning algorithm. It is called once on each node of the tree,
/// starting with the leaf nodes and working upwards. For each node, a vector of probabilities is added to the _probabilities
/// map. For leaf nodes, the probabilities are set based on input data, but for internal nodes the probabilities are calculated 
/// based on the probabilities of the child nodes. After all calculations, the caller should call the max_likelihood method,
/// passing the root of the tree, to determine the overall probabilities for the tree.
void compute_node_probability(const clade *node, const gene_family&_gene_family, const error_model*_p_error_model,
    std::map<const clade *, std::vector<double> >& _probabilities,
    int _max_root_family_size,
    int _max_parsed_family_size,
    const lambda* _lambda,
    const matrix_cache& _calc);

/* START: Uniform distribution */
std::vector<int> uniform_dist(int n_draws, int min, int max);
/* END: Uniform distribution - */

clademap<int> create_family(pvalue_parameters p, int root_family_size);
void set_weighted_random_family_size(const clade *node, clademap<int> *sizemap, const lambda *p_lambda, error_model *p_error_model, int max_family_size, const matrix_cache& cache);
std::vector<double> get_random_probabilities(pvalue_parameters p, int number_of_simulations, int root_family_size);
size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

double pvalue(double v, const vector<double>& conddist);

//! computes a pvalue for each family. Returns a vector of pvalues matching the list of families
std::vector<double> compute_pvalues(pvalue_parameters p, const std::vector<gene_family>& families, int number_of_simulations);
std::vector<double> compute_family_probabilities(pvalue_parameters p, const vector<clademap<int>>& sizes);

#endif
