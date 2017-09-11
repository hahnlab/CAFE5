#ifndef CORE_H
#define CORE_H

#include "clade.h"
#include "probability.h"

class clade; // core.cpp includes clade.h before including core.h

class process;


class gamma_bundle
{
	std::vector<process *> processes;
public:
	void add(process *p) {
		processes.push_back(p);
	}

	void prune();
};


class core {
private:
    std::ostream & _ost; 
    lambda *_p_lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _max_root_family_size;
    int _max_family_size_sim;
    int _total_n_families_sim;
    vector<gene_family> * _p_gene_families;
    vector<int> _rootdist_vec; // in case the user wants to use a specific root size distribution for all simulations
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin
    
    //! Gamma
    vector<double> _lambda_multipliers;
    vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category
    vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families
    double _alpha;
    
    //! Simulations
    vector<process*> _sim_processes; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
    
    //! Inference
	vector<gamma_bundle> _inference_bundles; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
public:
    core(): _ost(cout), _p_lambda(NULL), _p_tree(NULL), _p_gene_families(NULL), _total_n_families_sim(1), _lambda_multipliers(1) {}
    
    core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
		int n_gamma_cats, double alpha);
    
    core(ostream & ost, lambda* p_lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
		vector<int>& cats, vector<double>&mul);
    
    //void estimate_processes(); 
    
    //! Setter methods
    void set_lambda(lambda *p_lambda);
    
    void set_tree(clade *p_tree);
    
    void set_gene_families(std::vector<gene_family> *p_gene_families);
    
    void set_max_sizes(int max_family_size, int max_root_family_size);
    
    void set_max_size_sim(int max_family_size_sim);
    
    void set_rootdist_vec(std::vector<int> rootdist_vec);
    
    void set_total_n_families_sim(int total_n_families_sim);
    
    void set_alpha(double alpha);
    
    void set_lambda_multipliers(std::vector<double> lambda_multipliers);
    
    void set_lambda_bins(std::vector<int> lambda_bins);
    
    //! Simulation methods
    void start_sim_processes();

    void start_inference_processes();

    void simulate_processes();

    //! Inference methods
    void infer_processes();

    //! Gamma methods
    void adjust_n_gamma_cats(int n_gamma_cats);
    
    //! Printing methods
    void print_parameter_values();
    
    void print_simulations(ostream& ost);
};

#endif /* CORE_H */

