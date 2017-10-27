#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>

class clade;
class probability_calculator;
class gene_family;
class root_equilibrium_distribution;
class core;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

//! Abstract class to hold lambda value.
/*!
  The main feature of this abstract class is its pure virtual member calculate_child_factor. This method is required because one or more lambda may be specified, and this determines how the pruner should compute the likelihoods. So this abstract class is what decides how the likelihood is computed for each descendant.
 */
class lambda {
protected:
	probability_calculator *_p_calc;
public:
	lambda(probability_calculator *p_calc) : _p_calc(p_calc)
	{
			
	}
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) = 0; //!< Pure virtual function (= 0 is the 'pure specifier' and indicates this function MUST be overridden by a derived class' method)
    virtual lambda *multiply(double factor) = 0;
};

//! (lambda) Derived class 1: one lambda for whole tree
class single_lambda : public lambda {
private:
    double _lambda;
    
public:
    single_lambda(probability_calculator *p_calc, double lam) : lambda(p_calc), _lambda(lam) { } //!< Constructor 
    double get_single_lambda() const { return _lambda; }
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size); //!< Computes tr. prob. matrix, and multiplies by likelihood vector. Returns result (=factor).

	virtual lambda *multiply(double factor)
	{
		return new single_lambda(_p_calc, _lambda * factor);
	}
};

//! (lambda) Derived class 2: multiple lambdas
class multiple_lambda : public lambda {
private:
    std::map<std::string, int> _node_name_to_lambda_index;
    std::vector<double> _lambdas;
    
public:
    multiple_lambda(probability_calculator *p_calc, std::map<std::string, int> nodename_index_map, std::vector<double> lambda_vector) : lambda(p_calc),
		_node_name_to_lambda_index(nodename_index_map), _lambdas(lambda_vector) { } //!< Constructor
    virtual std::vector<double> calculate_child_factor(clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size); //!< Computes tr. prob. matrix (uses right lambda for each branch) and multiplies by likelihood vector. Returns result (=factor).
    virtual lambda *multiply(double factor)
    {
        auto npi = _lambdas;

        for (auto& i : npi)
            i *= factor;

        return new multiple_lambda(_p_calc, _node_name_to_lambda_index, npi);
    }
};

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

std::vector<double> get_posterior(std::vector<gene_family> gene_families, int max_family_size, int max_root_family_size, double lambda, clade *p_tree);
double find_best_lambda(core *p_model, root_equilibrium_distribution *p_distribution);

#endif
