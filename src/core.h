#ifndef CORE_H
#define CORE_H

#include "clade.h"

class clade; // core.cpp includes clade.h before including core.h

class process;

class core {
private:
    std::ostream & _ost; 
    double _lambda; // TODO: multiple lambdas for different branches
    clade *_p_tree;
    int _max_family_size;
    int _total_n_families;
    vector<int> _rootdist_vec;
    vector<double> _lambda_multipliers;
    vector<int> _lambda_bins; // must be same size as root family distribution total count
    vector<process*> _processes; // as of now, each process will be ONE simulation (i.e., simulate ONE gene family) under ONE lambda multiplier; maybe later each process will run more simulations under the same lambda multiplier?
    
public:
    core(): _ost(cout), _total_n_families(1), _lambda_multipliers(1) {}
    core(ostream & ost, double lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec, vector<double> lambda_multipliers, vector<int> lambda_bins): _ost(ost), _lambda(lambda), _p_tree(p_tree), _max_family_size(max_family_size), _total_n_families(total_n_families), _rootdist_vec(rootdist_vec), _lambda_multipliers(lambda_multipliers), _lambda_bins(lambda_bins) {} // Ben: I need to make this general enough so that if there is a single lambda multiplier, this also works
    void start_processes();
    void simulate_processes();
    //void estimate_processes(); 
    void print_parameter_values();
    void print_simulations(ostream& ost);
};

#endif /* CORE_H */

