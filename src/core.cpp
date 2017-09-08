#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <random>
#include "clade.h"
#include "core.h"
#include "family_generator.h"
#include "gamma.h"

/* START: Drawing random root size from uniform */
template<typename itr, typename random_generator>
itr select_randomly(itr start, itr end, random_generator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g)); // advances iterator (start) by dis(g), where g is a seeded mt generator
    return start;
}

template<typename itr>
itr select_randomly(itr start, itr end) {
    static std::random_device rd; // randomly generates uniformly distributed ints (seed)
    static std::mt19937 gen(rd()); // seeding Mersenne Twister generator
    return select_randomly(start, end, gen); // plug in mt generator to advance our container and draw random element from it
}
/* END: Drawing random root size from uniform */

class process {
private:
    ostream & _ost;
    lambda* _lambda;
    double _lambda_multiplier;
    clade *_p_tree;
    int _max_family_size;
    vector<int> _rootdist;
    int _root_size; // will be drawn from _rootdist by process itself
    trial *_my_simulation;
    
public:
    process(): _ost(cout), _lambda(NULL), _lambda_multiplier(1.0) {}
    
    process(ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree, int max_family_size, vector<int> rootdist): _ost(ost), _lambda(lambda), _lambda_multiplier(lambda_multiplier), _p_tree(p_tree), _max_family_size(max_family_size), _rootdist(rootdist) {

		if (_rootdist.empty())
		{
			// generating uniform root distribution when no distribution is provided 
			_rootdist.resize(max_family_size);
			for (size_t i = 0; i < _rootdist.size(); ++i)
				_rootdist[i] = i;
		}

		_root_size = *select_randomly(_rootdist.begin(), _rootdist.end()); // getting a random root size from the provided (core's) root distribution
        cout << "_root_size is " << _root_size << endl;
    }
    
    void run_simulation();
    
    void print_simulation(std::ostream & ost);
    
    trial * get_simulation();
};

//! Run process' simulation
void process::run_simulation() {
    double lambda_m = dynamic_cast<single_lambda*>(_lambda)->get_single_lambda() * _lambda_multiplier;
    _my_simulation = simulate_family_from_root_size(_p_tree, _root_size, _max_family_size, lambda_m);
}

//! Printing process' simulation
void process::print_simulation(std::ostream & ost) {

    // Printing gene counts
    for (trial::iterator it = _my_simulation->begin(); it != _my_simulation->end(); ++it) {
	ost << it->second << "\t";
    }

    ost << _lambda_multiplier << endl;
} 

//! Return simulation
trial * process::get_simulation() {
    return _my_simulation;
}

core::core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
	int n_gamma_cats, double alpha) : _ost(ost), _lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size), 
	_total_n_families(total_n_families), _rootdist_vec(rootdist_vec), _gamma_cat_probs(n_gamma_cats),
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
	vector<int>& cats, vector<double>&mul) : _ost(ost), _lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size),
	_total_n_families(total_n_families), _rootdist_vec(rootdist_vec),
	_gamma_cats(cats), _lambda_multipliers(mul)
{
}

//! Populate _processes (vector of processes)
void core::start_sim_processes() {
    
    for (int i = 0; i < _total_n_families; ++i) {
        double lambda_bin = _gamma_cats[i];      
        process *p_new_process = new process(_ost, _lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
        _sim_processes.push_back(p_new_process);
    }
    
    cout << _sim_processes.size() << " processes have been started." << endl;
}

//! Run simulations in all processes, in series... (TODO: in parallel!)
void core::simulate_processes() {
    for (int i  = 0; i < _total_n_families; ++i) {
        _sim_processes[i]->run_simulation();
    }
}

//! Print processes' simulations
void core::print_simulations(std::ostream& ost) {
    
    // Printing header
    for (trial::iterator it = _sim_processes[0]->get_simulation()->begin(); it != _sim_processes[0]->get_simulation()->end(); ++it) {
	ost << "#" << it->first->get_taxon_name() << endl;
    }
    
    for (int i = 0; i < _total_n_families; ++i) {
        _sim_processes[i]->print_simulation(cout);
    }
}

/* TODO: later this will become a member of the core class, which is the wrapper of the process class*/
void core::print_parameter_values() {
    
    cout << endl << "You have set the following parameter values:" << endl;
    
    if (dynamic_cast<single_lambda*>(_lambda)->get_single_lambda() == 0.0) {
        cout << "Lambda has not been specified." << endl;
    }
    else {
        cout << "Lambda: " << _lambda << endl;
    }
    
    if (_p_tree == NULL) {
        cout << "A tree has not been specified." << endl;
    }
    else {
        cout << "The tree is:" << endl;
        _p_tree->print_clade();
    }
    
}
