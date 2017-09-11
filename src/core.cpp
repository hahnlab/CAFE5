#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include "clade.h"
#include "core.h"
#include "family_generator.h"
#include "gamma.h"
#include "process.h"

void gamma_bundle::prune() {
    for (int i = 0; i < processes.size(); ++i)
	processes[i]->prune();
}

//! Simulation: core constructor when just alpha is provided.
core::core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
	int n_gamma_cats, double alpha) : _ost(ost), _p_lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size), 
	_total_n_families_sim(total_n_families), _rootdist_vec(rootdist_vec), _gamma_cat_probs(n_gamma_cats),
	_lambda_multipliers(n_gamma_cats)
{
	if (!rootdist_vec.empty()) {
		_rootdist_bins.push_back(rootdist_vec); // just 1st element
	}

	else {
		_rootdist_vec = uniform_dist(total_n_families, 1, max_family_size); // the user did not specify one... using uniform from 1 to max_family_size!
		_rootdist_bins.push_back(_rootdist_vec); // just 1st element (but we could specify different root dists for each lambda bin)
	}

	if (n_gamma_cats > 1)
	{
		get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference
		auto cats = weighted_cat_draw(total_n_families, _gamma_cat_probs);
		_gamma_cats = *cats;
		delete cats;
	}

	for (auto i = _gamma_cat_probs.begin(); i != _gamma_cat_probs.end(); ++i) {
		cout << "Should be all the same probability: " << *i << endl;
	}

	for (auto i = _lambda_multipliers.begin(); i != _lambda_multipliers.end(); ++i) {
		cout << "Lambda multiplier (rho) is: " << *i << endl;
	}

	for (auto i = _gamma_cats.begin(); i != _gamma_cats.end(); ++i) {
		cout << "Gamma category is: " << *i << endl;
	}

}

core::core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
	vector<int>& cats, vector<double>&mul) : _ost(ost), _p_lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size),
	_total_n_families_sim(total_n_families), _rootdist_vec(rootdist_vec),
	_gamma_cats(cats), _lambda_multipliers(mul)
{
}

//! Set pointer to lambda in core class
void core::set_lambda(lambda *p_lambda) {
    _p_lambda = p_lambda;
}

//! Set pointer to lambda in core class
void core::set_tree(clade *p_tree) {
    _p_tree = p_tree;
}

//! Set pointer to vector of gene family class instances
void core::set_gene_families(std::vector<gene_family> *p_gene_families) {
    _p_gene_families = p_gene_families;
}

//! Set max family sizes and max root family sizes for INFERENCE
void core::set_max_sizes(int max_family_size, int max_root_family_size) {
    _max_family_size = max_family_size;
    _max_root_family_size = max_root_family_size;
}

//! Set root distribution vector
void core::set_rootdist_vec(std::vector<int> rootdist_vec) {
    _rootdist_vec = rootdist_vec;
}

//! Set total number of families to simulate
void core::set_total_n_families_sim(int total_n_families_sim) {
    _total_n_families_sim = total_n_families_sim;
}

//! Resize all gamma-related vectors according to provided number (integer) of gamma categories
void core::adjust_n_gamma_cats(int n_gamma_cats) {
    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
}

//! Resize gamma_cats vector that assigns gamma class membership of families to be inferred/simulated
void core::adjust_family_gamma_membership(int n_families) {
    _gamma_cats.resize(n_families);
}

//! Set alpha for gamma distribution
void core::set_alpha(double alpha) {
    _alpha = alpha;
    if (_gamma_cats.size() > 1)
	get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference
}

//! Set lambda multipliers for each gamma category
void core::set_lambda_multipliers(std::vector<double> lambda_multipliers) {
    _lambda_multipliers = lambda_multipliers;
}

//! Set lambda bins (each int is a vector pointing to a gamma category)
void core::set_lambda_bins(std::vector<int> lambda_bins) {
    _gamma_cats = lambda_bins;
}

//! Populate _processes (vector of processes)
void core::start_sim_processes() {
    
    for (int i = 0; i < _total_n_families_sim; ++i) {
        double lambda_bin = _gamma_cats[i];      
        simulation_process *p_new_process = new simulation_process(_ost, _p_lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
        _sim_processes.push_back(p_new_process);
    }
    
    // cout << _sim_processes.size() << " processes have been started." << endl;
}

void gamma_core::start_inference_processes() {

    for (int i = 0; i < _p_gene_families->size(); ++i) {
	gamma_bundle bundle;
        
        cout << "Started inference bundle " << i+1 << endl;
	
        for (int j = 0; j < _gamma_cat_probs.size(); ++j) {
            double lambda_bin = _gamma_cat_probs[j];
			inference_process *p_new_process = new inference_process(_ost, _p_lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
            bundle.add(p_new_process);
            
            cout << "  Started inference process " << j+1 << endl;
            
	}

	_inference_bundles.push_back(bundle);
    }
}


//! Run simulations in all processes, in series... (TODO: in parallel!)
void core::simulate_processes() {
    for (int i  = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->run_simulation();
    }
}

void core::start_inference_processes()
{
}

void core::infer_processes()
{
	// prune all the families with the same lambda
}

//! Infer bundle
void gamma_core::infer_processes() {
    for (int i = 0; i < _inference_bundles.size(); ++i) {
        cout << endl << "About to prune a gamma bundle." << endl;
	_inference_bundles[i].prune();
    }
}

//! Print processes' simulations
void core::adjust_family(std::ostream& ost) {
    
    // Printing header
    for (trial::iterator it = _sim_processes[0]->get_simulation()->begin(); it != _sim_processes[0]->get_simulation()->end(); ++it) {
	ost << "#" << it->first->get_taxon_name() << endl;
    }
    
    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->print_simulation(cout);
    }
}

/* TODO: later this will become a member of the core class, which is the wrapper of the process class*/
void core::print_parameter_values() {
    
    cout << endl << "You have set the following parameter values:" << endl;
    
    if (dynamic_cast<single_lambda*>(_p_lambda)->get_single_lambda() == 0.0) {
        cout << "Lambda has not been specified." << endl;
    }
    else {
        cout << "Lambda: " << _p_lambda << endl;
    }
    
    if (_p_tree == NULL) {
        cout << "A tree has not been specified." << endl;
    }
    else {
        cout << "The tree is:" << endl;
        _p_tree->print_clade();
    }
    
}
