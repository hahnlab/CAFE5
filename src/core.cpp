#include <vector>
#include <iostream>
#include "clade.h"
#include "core.h"

class process {
private:
    double _lambda;
    double _lambda_multiplier;
    clade *_p_tree;
    int _max_family_size;
    int _n_simulations;
    
public:
    process(): _lambda(0.0), _lambda_multiplier(1.0) {}  
    process(double lambda, double lambda_multiplier, clade *p_tree, int max_family_size) {}
    void run_simulations();
};

//! Populate _processes (vector of processes)
void core::start_processes() {
    for (int i = 0; i < _n_processes; ++i) {
        process *p_new_process = new process(_lambda, _lambda_multipliers[i], _p_tree, _max_family_size); // if a single _lambda_multiplier, how do we do it?
        _processes.push_back(p_new_process);
    }
}

/* TODO: later this will become a member of the core class, which is the wrapper of the process class*/
void core::print_parameter_values() {
    
    cout << endl << "You have set the following parameter values:" << endl;
    
    if (_lambda == 0.0) {
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
    
    cout << "There is/are " << _n_processes << " rate class(es)." << endl;
    
}