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
#include "poisson.h"

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
    return new simulation_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec, family_number); // if a single _lambda_multiplier, how do we do it?
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

base_core::~base_core()
{
    for (size_t i = 0; i < processes.size(); ++i)
        delete processes[i];
}

void base_core::start_inference_processes()
{
    processes.clear();
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

double base_core::infer_processes(prior_distribution *prior) {
    initialize_rootdist_if_necessary();
    prior->initialize(_rootdist_vec);

    results.clear();
    std::vector<double> all_families_likelihood(processes.size());
    // prune all the families with the same lambda
    for (int i = 0; i < processes.size(); ++i) {
        std::cout << "Process " << i << std::endl;
        auto partial_likelihood = processes[i]->prune();
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = prior->compute(j);
//            std::cout << "log-eq_prob = " << std::log(eq_freq) << ", partial log-lk = " << std::log(partial_likelihood[j]) << std::endl;

            double log_full_lk = std::log(partial_likelihood[j]) + std::log(eq_freq);
            full[j] = log_full_lk;

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }
        
        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)

        results.push_back(family_info_stash(i, 0.0, 0.0, 0.0, all_families_likelihood[i], false));

        std::cout << "lnL of family " << i << ": " << all_families_likelihood[i] << std::endl;
    }

    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    std::cout << "-lnL: " << final_likelihood << std::endl;

    return final_likelihood;
}

void base_core::print_results(std::ostream& ost)
{
    ost << "#FamilyID\tLikelihood of Family" << endl;
    for (const auto& r : results)
    {
        ost << r.family_id << "\t" << r.posterior_probability << endl;
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

class max_branch_length_finder
{
    double _result;
public:
    max_branch_length_finder() : _result(0.0)
    {

    }
    void operator()(clade *c)
    {
        if (c->get_branch_length() > _result)
            _result = c->get_branch_length();
    }
    double result() const {
        return _result;
    }
};


double core::initialize_lambda_guess()
{
    max_branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);
    return finder.result();
}

float equilibrium_frequency::compute(int val) const
{
    int sum = std::accumulate(_rootdist_vec.begin(), _rootdist_vec.end(), 0);
    return float(_rootdist_vec[val]) / float(sum);
}

void poisson_frequency::initialize(std::vector<int> rootdist_vec)
{
    poisson = get_prior_rfsize_poisson_lambda(0, rootdist_vec.size(), _poisson_lambda);
}
