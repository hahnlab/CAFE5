#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <iostream>
#include "utils.h"
#include "io.h"

class clade;
double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size);
double chooseln(double n, double k);
double unifrnd();

class lambda
{
public:
	virtual std::vector<double> calculate_child_factor(clade *child, std::size_t sz, std::vector<double> probabilities) = 0;
	// = 0 implies that this class does NOT implement this method, but any descendant class MUST implement it
};

class single_lambda : public lambda
{
	double _lambda;
public:
	single_lambda(double lambda) : _lambda(lambda)
	{

	}

	double get_single_lambda() const { return _lambda; }
	virtual std::vector<double> calculate_child_factor(clade *child, std::size_t sz, std::vector<double> probabilities);
};

class multiple_lambda : public lambda
{
	std::map<std::string, int> _node_name_to_lambda_index;
	std::vector<double> _lambdas;
public:
	multiple_lambda(std::map<std::string, int> a, std::vector<double> b) : _node_name_to_lambda_index(a), _lambdas(b)
	{
	}
	virtual std::vector<double> calculate_child_factor(clade *child, std::size_t sz, std::vector<double> probabilities);
};


/* START: Likelihood computation ---------------------- */

class likelihood_computer {
private:
    std::map<clade *, std::vector<double> > _probabilities; //!< represents probability of the node having various family sizes
    gene_family *_family;
    int _max_possible_family_size;
    lambda* _lambda;
    
public:
    likelihood_computer(int max_possible_family_size, lambda* lambda, gene_family *family) : _max_possible_family_size(max_possible_family_size), _lambda(lambda) {
        _family = family;
    }
  
    void operator()(clade *node);

    std::vector<double> get_likelihoods(clade *node) const { 
        return _probabilities.at(node);
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
};

/* END: Probability calculator for simulator -------- */

/* START: Uniform distribution */
std::vector<int> uniform_dist(int n_draws, int min, int max);
/* END: Uniform distribution - */

std::vector<int> * weighted_cat_draw(int n_draws, std::vector<double> gamma_cat_probs);

#endif