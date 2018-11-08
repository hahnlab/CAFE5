#ifndef PROCESS_H
#define PROCESS_H

#include <iosfwd>
#include <vector>
#include <map>

class lambda;
class clade;
class gene_family;
class matrix_cache;
class root_equilibrium_distribution;
class child_multiplier;
class error_model;
class root_distribution;

typedef std::map<const clade *, int> trial;

class process {
protected:
	std::ostream & _ost;
	lambda* _lambda;
	double _lambda_multiplier;
	const clade *_p_tree;
	int _max_root_family_size;

	process(std::ostream & ost, lambda* lambda, double lambda_multiplier, const clade *p_tree, int max_root_family_size) :
		_ost(ost), _lambda(lambda), _lambda_multiplier(lambda_multiplier), _p_tree(p_tree),
		_max_root_family_size(max_root_family_size)
		{}
public:
    double get_lambda_multiplier() const {
        return _lambda_multiplier;
    }
};

/// Called when we fix lambda to calculate the probabilities of the tree
class inference_process : public process {
	const gene_family *_p_gene_family;
    const error_model *_p_error_model;
    int _max_family_size;
public:
	inference_process(std::ostream & ost, 
        lambda* lambda, 
        double lambda_multiplier, 
        const clade *p_tree, 
        int max_family_size,
		int max_root_family_size, 
        const gene_family *fam, 
        error_model *p_error_model) : 
        process(ost, lambda, lambda_multiplier, p_tree,
			max_root_family_size),
            _max_family_size(max_family_size),
            _p_error_model(p_error_model),
            _p_gene_family(fam) 
    {
	}

	std::vector<double> prune(matrix_cache& calc);    // returns likelihood of the tree for each family size
    std::string family_id() const;
};

class simulation_process : public process {
    error_model *_p_error_model;
	int _max_family_size_sim;
    int _root_size; // will be drawn from _rootdist_vec by process itself
public:
	simulation_process(std::ostream &ost, lambda* lambda, double lambda_multiplier, const clade *p_tree,
        int max_root_family_size, int max_family_size_sim, int root_size, error_model* p_error_model);

	trial* run_simulation(const matrix_cache& cache);
};

#endif
