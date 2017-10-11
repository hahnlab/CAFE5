#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include "clade.h"
#include "core.h"
#include "family_generator.h"
#include "process.h"

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r)
{
    ost << r.family_id << "\t" << r.lambda_multiplier << "\t" << r.category_likelihood << "\t" << r.family_likelihood;
    ost << "\t" << r.posterior_probability << "\t" << (r.significant ? "*" : "N/S");
    return ost;
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

void core::start_sim_processes() {

    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes.push_back(create_simulation_process(i));
    }

    // cout << _sim_processes.size() << " processes have been started." << endl;
}

simulation_process* base_core::create_simulation_process(int family_number) {
    return new simulation_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
}



//! Run simulations in all processes, in series... (TODO: in parallel!)
void core::simulate_processes() {
    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->run_simulation();
    }
}

void core::print_processes(std::ostream& ost) {

    trial *sim = _sim_processes[0]->get_simulation();
    for (trial::iterator it = sim->begin(); it != sim->end(); ++it) {
          ost << "#" << it->first->get_taxon_name() << "\n";
    }

    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->print_simulation(ost);
    }
}

float equilibrium_frequency::compute(int val) const
{
    int sum = std::accumulate(_rootdist_vec.begin(), _rootdist_vec.end(), 0);
    return float(_rootdist_vec[val]) / float(sum);
}

base_core::~base_core()
{
    for (size_t i = 0; i < processes.size(); ++i)
        delete processes[i];
}

void base_core::start_inference_processes()
{
    for (int i = 0; i < _p_gene_families->size(); ++i) {

        cout << "Started inference process " << i + 1 << endl;

            //            double lambda_bin = _gamma_cat_probs[j];
        inference_process *p_new_process = new inference_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
        processes.push_back(p_new_process);
    }
}

void core::initialize_rootdist_if_necessary()
{
    if (_rootdist_vec.empty())
    {
        _rootdist_vec.resize(_max_root_family_size);
        std::fill(_rootdist_vec.begin(), _rootdist_vec.end(), 1);
    }

}

double base_core::infer_processes() {
    initialize_rootdist_if_necessary();
    equilibrium_frequency eq(_rootdist_vec);
    std::vector<double> all_families_likelihood(processes.size());
    // prune all the families with the same lambda
    for (int i = 0; i < processes.size(); ++i) {
        cout << "Process " << i << endl;
        auto partial_likelihood = processes[i]->prune();
        std::vector<double> full(partial_likelihood.size());
        
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = eq.compute(j);
            cout << "log-eq_prob = " << std::log(eq_freq) << ", partial log-lk = " << std::log(partial_likelihood[j]) << endl;
            
            double log_full_lk = std::log(partial_likelihood[j]) + std::log(eq_freq);
//            full[j] = log_full_lk;
            if (!isinf(log_full_lk))
                full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
            else
                full[j] = 0.0;         
            // cout << "Full lk of size " << j << ": " << full[j] << endl;
        }
        
        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
//        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
        cout << "lnL of family " << i << ": " << all_families_likelihood[i] << endl;
    }

    double multi = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    cout << "-lnL: " << multi << std::endl;

    return multi;
}

void base_core::print_results(std::ostream& ost)
{
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
