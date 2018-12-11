#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>

class clade;
class matrix_cache;
class gene_family;
class root_equilibrium_distribution;
class model;

struct FMinSearch;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

//! Abstract class to hold lambda value.
/*!
  The main feature of this abstract class is its pure virtual member calculate_child_factor. This method is required because one or more lambda may be specified, and this determines how the pruner should compute the likelihoods. So this abstract class is what decides how the likelihood is computed for each descendant.
 */
class lambda {
public:
    virtual std::vector<double> calculate_child_factor(const matrix_cache& calc, const clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const = 0; //!< Pure virtual function (= 0 is the 'pure specifier' and indicates this function MUST be overridden by a derived class' method)
    virtual lambda *multiply(double factor) const = 0;
    virtual void update(double* values) = 0;
    virtual int count() const = 0;
    virtual std::string to_string() const = 0;
    virtual double get_value_for_clade(const clade *c) const = 0;
    virtual bool is_valid() = 0;
    virtual lambda* clone() const = 0;

    virtual ~lambda() {}
};

//! (lambda) Derived class 1: one lambda for whole tree
class single_lambda : public lambda {
private:
    double _lambda;
    
public:
    single_lambda(double lam) : _lambda(lam) { } //!< Constructor 
    double get_single_lambda() const { return _lambda; }
    virtual std::vector<double> calculate_child_factor(const matrix_cache& calc, const clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const; //!< Computes tr. prob. matrix, and multiplies by likelihood vector. Returns result (=factor).

	virtual lambda *multiply(double factor) const
	{
		return new single_lambda(_lambda * factor);
	}
    virtual void update(double* values) { _lambda = *values; }

    virtual int count() const {
        return 1;
    }
    virtual std::string to_string() const;
    virtual double get_value_for_clade(const clade *c) const {
        return _lambda;
    }
    virtual bool is_valid() {
        return _lambda > 0;
    }
    virtual lambda *clone() const {
        return new single_lambda(_lambda);
    }
};

//! (lambda) Derived class 2: multiple lambdas
class multiple_lambda : public lambda {
private:
    std::map<std::string, int> _node_name_to_lambda_index;
    std::vector<double> _lambdas;
    
public:
    multiple_lambda(std::map<std::string, int> nodename_index_map, std::vector<double> lambda_vector) :
		_node_name_to_lambda_index(nodename_index_map), _lambdas(lambda_vector) { } //!< Constructor
    virtual std::vector<double> calculate_child_factor(const matrix_cache& calc, const clade *child, std::vector<double> probabilities, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const; //!< Computes tr. prob. matrix (uses right lambda for each branch) and multiplies by likelihood vector. Returns result (=factor).
    virtual lambda *multiply(double factor) const
    {
        auto npi = _lambdas;

        for (auto& i : npi)
            i *= factor;

        return new multiple_lambda(_node_name_to_lambda_index, npi);
    }
    virtual void update(double* values);    
    virtual int count() const {
        return _lambdas.size();
    }
    virtual std::string to_string() const;
    virtual double get_value_for_clade(const clade *c) const;
    virtual bool is_valid();

    std::vector<double> get_lambdas() const {
        return _lambdas;
    }
    virtual lambda *clone() const {
        return new multiple_lambda(_node_name_to_lambda_index, _lambdas);
    }

};

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

class optimizer_scorer;

class optimizer {
    FMinSearch* pfm;
    optimizer_scorer *_p_scorer;
public:
    optimizer(optimizer_scorer *scorer);

    struct result {
        std::vector<double> values;
        double score;
        int num_iterations;
    };

    result optimize();

    void log_results(FMinSearch * pfm, std::vector<double> &initial, double * re);

    bool quiet = false;
    bool explode = false;

};

inline std::ostream& operator<<(std::ostream& ost, const lambda& lambda)
{
    ost << lambda.to_string();
    return ost;
}

#endif
