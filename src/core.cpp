#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include "clade.h"
#include "core.h"
#include "family_generator.h"

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

        _root_size = *select_randomly(_rootdist.begin(), _rootdist.end()); // getting a random root size from the provided (core's) root distribution
        //cout << "_root_size is " << _root_size << endl;
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

//! Populate _processes (vector of processes)
void core::start_processes() {
    
    for (int i = 0; i < _total_n_families; ++i) {
        double lambda_bin = _lambda_bins[i];      
        process *p_new_process = new process(_ost, _lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
        _processes.push_back(p_new_process);
    }
    
    cout << _processes.size() << " processes have been started." << endl;
}

//! Run simulations in all processes, in series... (TODO: in parallel!)
void core::simulate_processes() {
    for (int i  = 0; i < _total_n_families; ++i) {
        _processes[i]->run_simulation();
    }
}

//! Print processes' simulations
void core::print_simulations(std::ostream& ost) {
    
    // Printing header
    for (trial::iterator it = _processes[0]->get_simulation()->begin(); it != _processes[0]->get_simulation()->end(); ++it) {
	ost << "#" << it->first->get_taxon_name() << endl;
    }
    
    for (int i = 0; i < _total_n_families; ++i) {
        _processes[i]->print_simulation(cout);
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
