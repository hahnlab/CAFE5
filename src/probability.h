#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <tuple>
#include <set>

#include <iostream>
#include <set>
#include "utils.h"
#include "io.h"
#include "lambda.h"

class clade;
class matrix;
class matrix_cache;

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff);
double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size);
double chooseln(double n, double k);
double unifrnd();

/* START: Likelihood computation ---------------------- */

/// The likelihood_computer class is a functor that is used in the pruning algorithm. It is called once on each node of the tree,
/// starting with the leaf nodes and working upwards. For each node, a vector of probabilities is added to the _probabilities
/// map. For leaf nodes, the probabilities are set based on input data, but for internal nodes the probabilities are calculated 
/// based on the probabilities of the child nodes. After all calculations, the caller should call the max_likelihood method,
/// passing the root of the tree, to determine the overall probabilities for the tree.
class likelihood_computer {
private:
    std::map<const clade *, std::vector<double> > _probabilities; //!< represents probability of the node having various family sizes
    int _max_root_family_size;
	int _max_parsed_family_size;
	const lambda* _lambda;
    const matrix_cache& _calc;
    const gene_family& _gene_family;
    const error_model* _p_error_model;
    
public:
    likelihood_computer(int max_root_family_size, int max_parsed_family_size, const lambda* lambda, const gene_family& family,
        const matrix_cache& calc, const error_model *p_error_model) : 
		_max_root_family_size(max_root_family_size),
		_max_parsed_family_size(max_parsed_family_size),
        _gene_family(family),
		_lambda(lambda),
        _calc(calc),
        _p_error_model(p_error_model) {
    }
  
    void operator()(const clade *node);

    std::vector<double> get_likelihoods(const clade *node) const { 
        return _probabilities.at(node);
    }

	double max_likelihood(const clade *node) const
	{
		if (_probabilities.at(node).empty())
			throw std::runtime_error("No probabilities calculated");

		// use "at" rather than [] so the method can be const
		return *std::max_element(_probabilities.at(node).begin(), _probabilities.at(node).end());
	}

};

/* END: Likelihood computation ---------------------- */

/* START: Uniform distribution */
std::vector<int> uniform_dist(int n_draws, int min, int max);
/* END: Uniform distribution - */

class category_selector
{
    std::vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families
public:
    /* START: Weighted draw from vector */
    //! Draw ints or doubles from n_draws equal intervals using specified weights.
    void weighted_cat_draw(int n_draws, std::vector<double> gamma_cat_probs);

    bool empty() const { return _gamma_cats.empty(); }

    int draw(int n) { return _gamma_cats[n]; }
};

std::vector<std::vector<double> > get_conditional_distribution_matrix(const clade *p_tree, int root_family_size, int max_family_size, int number_of_simulations, const lambda *p_lambda, const matrix_cache& cache);

#endif