#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>

class clade;

class gene_family;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

//! Abstract class to hold lambda value.
/*!
  The main feature of this abstract class is its pure virtual member calculate_child_factor. This method is required because one or more lambda may be specified, and this determines how the pruner should compute the likelihoods. So this abstract class is what decides how the likelihood is computed for each descendant.
 */
class lambda {
public:
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) = 0; //!< Pure virtual function (= 0 is the 'pure specifier' and indicates this function MUST be overridden by a derived class' method)
};

//! (lambda) Derived class 1: one lambda for whole tree
class single_lambda : public lambda {
private:
    double _lambda;
    
public:
    single_lambda(double lambda) : _lambda(lambda) { } //!< Constructor 
    double get_single_lambda() const { return _lambda; }
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size); //!< Computes tr. prob. matrix, and multiplies by likelihood vector. Returns result (=factor).
};

//! (lambda) Derived class 2: multiple lambdas
class multiple_lambda : public lambda {
private:
    std::map<std::string, int> _node_name_to_lambda_index;
    std::vector<double> _lambdas;
    
public:
    multiple_lambda(std::map<std::string, int> nodename_index_map, std::vector<double> lambda_vector) : _node_name_to_lambda_index(nodename_index_map), _lambdas(lambda_vector) { } //!< Constructor
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size); //!< Computes tr. prob. matrix (uses right lambda for each branch) and multiplies by likelihood vector. Returns result (=factor).
};

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

struct lambda_search_params
{
	lambda_search_params(clade *pt, std::vector<gene_family> f, int mfs, int rfms) : ptree(pt), families(f), 
		max_family_size(mfs), max_root_family_size(rfms), initial_lambda(0)
	{

	}
	clade *ptree;
	std::vector<gene_family> families;
	int max_family_size;
	int max_root_family_size;
	double initial_lambda;
};


std::vector<double> get_posterior(std::vector<gene_family> gene_families, int max_family_size, int max_root_family_size, double lambda, clade *p_tree);
double find_best_lambda(lambda_search_params *params);

#endif
