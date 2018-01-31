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

class likelihood_computer {
private:
    std::map<clade *, std::vector<double> > _probabilities; //!< represents probability of the node having various family sizes
    gene_family *_family;
    int _max_root_family_size;
	int _max_parsed_family_size;
	lambda* _lambda;
    matrix_cache& _calc;
    const error_model* _p_error_model;
    
public:
    likelihood_computer(int max_root_family_size, int max_parsed_family_size, lambda* lambda, gene_family *family,
        matrix_cache& calc, const error_model *p_error_model) : 
		_max_root_family_size(max_root_family_size),
		_max_parsed_family_size(max_parsed_family_size),
		_lambda(lambda),
        _calc(calc),
        _p_error_model(p_error_model) {
        _family = family;
    }
  
    void operator()(clade *node);

    std::vector<double> get_likelihoods(clade *node) const { 
        return _probabilities.at(node);
    }

	double max_likelihood(clade *node) const
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

std::vector<int> * weighted_cat_draw(int n_draws, std::vector<double> gamma_cat_probs);

std::vector<std::vector<double> > get_conditional_distribution_matrix(clade *p_tree, int root_family_size, int max_family_size, int number_of_simulations, double lambda);

#endif