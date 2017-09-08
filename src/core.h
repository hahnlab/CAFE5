#ifndef CORE_H
#define CORE_H

#include "clade.h"
#include "probability.h"

class clade; // core.cpp includes clade.h before including core.h

class process;

class core {
private:
    std::ostream & _ost; 
    lambda *_lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _total_n_families;
    vector<int> _rootdist_vec; // in case the user wants to use a specific root size distribution for all simulations
    vector<vector<int> > _rootdist_bins; // holds the distribution for each lambda bin
    vector<double> _lambda_multipliers;
	vector<double> _gamma_cat_probs;	// 
	vector<int> _gamma_cats; // each item is an index to a gamma category, from 0 to n_cat; vector must be of length = _total_n_families
    vector<process*> _sim_processes; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier
    
public:
    core(): _ost(cout), _total_n_families(1), _lambda_multipliers(1) {}
    
	core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
		int n_gamma_cats, double alpha);
    
	core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
		vector<int>& cats, vector<double>&mul);

    void start_sim_processes();
    
    void simulate_processes();
    
    //void estimate_processes(); 
    
    void print_parameter_values();
    
    void print_simulations(ostream& ost);
};

#endif /* CORE_H */

