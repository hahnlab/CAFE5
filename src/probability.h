#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <iostream>
#include "utils.h"
#include "io.h"
#include "lambda.h"

class clade;
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
    
public:
    likelihood_computer(int max_root_family_size, int max_parsed_family_size, lambda* lambda, gene_family *family) : 
		_max_root_family_size(max_root_family_size),
		_max_parsed_family_size(max_parsed_family_size),
		_lambda(lambda) {
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

/* START: Probability calculator for simulator ------ */

//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
  Contains a map (_cache) that serves as a hash table to store precalculated values.
  If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class probability_calculator {
private:
    // C++ idiomatic way of creating a key that can be compared to another key
    struct key {
        double lambda;
        double branch_length;
        int parent_size;
        int size;

        // std::tie : when given n input values, it returns a std:tuple (a container) of size n
        // the < operator is defined for std::tuple, and it allows lexicographical sorting of n-tuples (it sorts and returns the first); the std::tuple < operator is left intact below -- below we overload the < operator of the probability_calculator class
        //! Operator < overload
        /*!
          Necessary for the map method find() used below
        */ 
        bool operator<(const key &o) const {
            return std::tie(size, parent_size, branch_length, lambda) < std::tie(o.size, o.parent_size, o.branch_length, o.lambda);
        }
    };
  
    std::map<key, double> _cache; //!< map that stores transition probabilities (given a parent and child size, branch length and lambda)
  
public:
    double get_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size) {
        key k = { lambda, branch_length, parent_size, size };
    
        if (_cache.find(k) == _cache.end()) { // if k is not in _cache
            _cache[k] = the_probability_of_going_from_parent_fam_size_to_c(lambda, branch_length, parent_size, size);
        }

        return _cache[k];
    }
    
    void print_cache(std::ostream &ost) const {
        for (auto key_value = _cache.begin(); key_value != _cache.end(); ++key_value) {
            ost << key_value.first << ", " << key_value.second << endl;
        }
    }
};

vector<vector<double> > get_matrix(int size, int branch_length, double lambda);

vector<double> matrix_multiply(const vector<vector<double> >& matrix, const vector<double>& v, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size);

/* END: Probability calculator for simulator -------- */

/* START: Uniform distribution */
std::vector<int> uniform_dist(int n_draws, int min, int max);
/* END: Uniform distribution - */

std::vector<int> * weighted_cat_draw(int n_draws, std::vector<double> gamma_cat_probs);

std::vector<std::vector<double> > get_conditional_distribution_matrix(clade *p_tree, int root_family_size, int max_family_size, int number_of_simulations, double lambda);

#endif