#ifndef CORE_H
#define CORE_H

#include "clade.h"
#include "probability.h"

class clade; // core.cpp includes clade.h before including core.h

class process;
class inference_process;
class simulation_process;

class gamma_bundle {
	std::vector<inference_process *> processes;
public:
	void add(inference_process *p) {
		processes.push_back(p);
	}

	void prune();
};


class core {
protected:
    std::ostream & _ost; 
    lambda *_p_lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _max_root_family_size;
    int _total_n_families_sim;
    vector<gene_family> * _p_gene_families;
    vector<int> _rootdist_vec; // in case the user wants to use a specific root size distribution for all simulations
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin
    
    //! Simulations
    vector<simulation_process*> _sim_processes; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
    
public:
    //! Basic constructor
    core(): _ost(cout), _p_lambda(NULL), _p_tree(NULL), _p_gene_families(NULL), _total_n_families_sim(1) {}
    
    core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec) :         _ost(ost), _p_lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size), _total_n_families_sim(total_n_families), _rootdist_vec(rootdist_vec)
    {}
    
    //void estimate_processes(); 
    
    //! Setter methods
    void set_lambda(lambda *p_lambda);
    
    void set_tree(clade *p_tree);
    
    void set_gene_families(std::vector<gene_family> *p_gene_families);
    
    void set_max_sizes(int max_family_size, int max_root_family_size);
    
    void set_rootdist_vec(std::vector<int> rootdist_vec);
    
    void set_total_n_families_sim(int total_n_families_sim);
    
    //! Simulation methods
    void start_sim_processes();

    virtual simulation_process* create_simulation_process(int family_number) = 0;

    void simulate_processes();

    void print_processes(std::ostream& ost);

    //! Inference methods
    virtual void start_inference_processes() = 0;
    
    virtual void infer_processes() = 0;
    
    //! Printing methods
    void print_parameter_values();
    
    void adjust_family(ostream& ost);
};

class base_core : public core {
    std::vector<inference_process *> processes;
    virtual void start_inference_processes();
    virtual void infer_processes();
    virtual simulation_process* create_simulation_process(int family_number);
};

class gamma_core : public core {
private:
    //! Gamma
    vector<double> _lambda_multipliers;
    
    vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category
    
    vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families
    
    double _alpha;
    
    vector<gamma_bundle> _inference_bundles; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
    
public:
    //! Basic constructor
    gamma_core() { _alpha = 0; };
    
    //! Simulation constructor
    gamma_core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, int n_gamma_cats, double alpha);
    
    //! Simulation constructor
    gamma_core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, vector<int>& cats, vector<double>&mul);
    
    //! Gamma methods
    void adjust_n_gamma_cats(int n_gamma_cats);
    
    void adjust_family_gamma_membership(int n_families);
        
    //! Setters
    void set_alpha(double alpha, int n_families);
    void initialize_with_alpha(int n_gamma_cats, int n_families, double alpha);
    void initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats);

    void set_lambda_multipliers(std::vector<double> lambda_multipliers);
    
    void set_gamma_cats(std::vector<int> gamma_cats);
        
    //! Simulation methods
    virtual simulation_process* create_simulation_process(int family_number);

    //! Inference methods
    void start_inference_processes();

    void infer_processes();
};
#endif /* CORE_H */

