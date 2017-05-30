#ifndef CORE_H
#define CORE_H

class clade; // core.cpp includes clade.h before including core.h

class process;

class core {
private:
    double _lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _n_processes;
    //vector<int> _n_simulations_in_processes; // Ben: the sum of all elements in this vector is the total # of simulated families; each element tells core how many simulations each process should carry out (what process' _n_simulations should be)
    vector<double> _lambda_multipliers;
    vector<process*> _processes;
    
public:
    core(): _n_processes(1), _lambda_multipliers(1) {}
    core(double lambda, clade *p_tree, int max_family_size, int n_processes, vector<double> lambda_multipliers): _lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size), _n_processes(n_processes), _lambda_multipliers(lambda_multipliers) {} // Ben: I need to make this general enough so that if there is a single lambda multiplier, this also works
    void start_processes();
    void print_parameter_values();
    //void run_simulations(); // Ben: this function should run the specified number of simulations for each process, as specified in _n_simulations_in_processes
    
};

#endif /* CORE_H */

