#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include "clade.h"
#include "core.h"
#include "family_generator.h"
#include "gamma.h"
#include "process.h"

std::vector<family_info_stash> gamma_bundle::prune(const vector<double>& _gamma_cat_probs, equilibrium_frequency *eq) {
    assert(_gamma_cat_probs.size() == processes.size());

    std::vector<family_info_stash> result;
//    std::vector<double> all_gamma_cats_likelihood(processes.size());

    for (int k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        auto partial_likelihood = processes[k]->prune();
        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = eq->compute(j);
            full[j] = partial_likelihood[j] * eq_freq;
            // cout << "Likelihood " << j+1 << ": Partial " << partial_likelihood[j] << ", eq freq: " << eq_freq << ", Full " << full[j] << endl;
        }
        //all_gamma_cats_likelihood[k] = accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k];
        //cout << "Likelihood of gamma cat " << k << " = " << all_gamma_cats_likelihood[k] << std::endl;

        family_info_stash r;
        r.category_likelihood = accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k];
        r.lambda_multiplier = processes[k]->get_lambda_multiplier();
        result.push_back(r);
        
    }


    return result;
}

//! Simulation: gamma_core constructor when just alpha is provided.
gamma_core::gamma_core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
	int n_gamma_cats, double alpha) : core(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
        _gamma_cat_probs(n_gamma_cats), _lambda_multipliers(n_gamma_cats) {
	
    if (!rootdist_vec.empty()) {
        _rootdist_bins.push_back(rootdist_vec); // just 1st element
    }

    else {
        _rootdist_vec = uniform_dist(total_n_families, 1, max_family_size); // the user did not specify one... using uniform from 1 to max_family_size!
	_rootdist_bins.push_back(_rootdist_vec); // just 1st element (but we could specify different root dists for each lambda bin)
    }

    if (n_gamma_cats > 1) {
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

//! Simulation: core constructor when alpha is not provided and membership is provided directly, in addition to the lambda multipliers (something the user decides to use)
gamma_core::gamma_core(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
	vector<int>& cats, vector<double>&mul) : core(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
	_gamma_cats(cats), _lambda_multipliers(mul) {}

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

void gamma_core::initialize_with_alpha(int n_gamma_cats, int n_families, double alpha)
{
    adjust_n_gamma_cats(n_gamma_cats);
    adjust_family_gamma_membership(n_families);
    set_alpha(alpha, n_families);
}

void gamma_core::initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats)
{
    assert(lambda_multipliers.size() == n_gamma_cats);
    assert(gamma_cats.size() == n_families);
    //    adjust_n_gamma_cats(n_gamma_cats);
//    adjust_family_gamma_membership(n_families);
//    set_alpha(alpha, n_families);

    set_lambda_multipliers(lambda_multipliers);
    set_gamma_cats(gamma_cats);
}

//! Resize all gamma-related vectors according to provided number (integer) of gamma categories
void gamma_core::adjust_n_gamma_cats(int n_gamma_cats) {
    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
}

//! Resize gamma_cats vector that assigns gamma class membership of families to be inferred/simulated
void gamma_core::adjust_family_gamma_membership(int n_families) {
    _gamma_cats.resize(n_families);
}

//! Set alpha for gamma distribution
void gamma_core::set_alpha(double alpha, int n_families) {   
    _alpha = alpha;
    if (_gamma_cats.size() > 1)
    	get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference
    
    vector<int>* cats = weighted_cat_draw(n_families, _gamma_cat_probs);
    _gamma_cats = *cats;
    delete cats;

    for (std::vector<double>::iterator it = _gamma_cat_probs.begin(); it != _gamma_cat_probs.end(); ++it) {
        cout << "Gamma cat prob is : " << *it << endl;
    }
}

//! Set lambda multipliers for each gamma category
void gamma_core::set_lambda_multipliers(std::vector<double> lambda_multipliers) {
    _lambda_multipliers = lambda_multipliers;
}

//! Set lambda bins (each int is a vector pointing to a gamma category)
void gamma_core::set_gamma_cats(std::vector<int> gamma_cats) {
    _gamma_cats = gamma_cats;
}

void core::start_sim_processes() {

    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes.push_back(create_simulation_process(i));
    }

    // cout << _sim_processes.size() << " processes have been started." << endl;
}

//! Populate _processes (vector of processes)
simulation_process* gamma_core::create_simulation_process(int family_number) {
    double lambda_bin = _gamma_cats[family_number];
    return new simulation_process(_ost, _p_lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
}

simulation_process* base_core::create_simulation_process(int family_number) {
    return new simulation_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
}



void gamma_core::start_inference_processes() {

    for (int i = 0; i < _p_gene_families->size(); ++i) {
	    gamma_bundle bundle;
        
        cout << "Started inference bundle " << i+1 << endl;
	
        for (int j = 0; j < _lambda_multipliers.size(); ++j) {
//            double lambda_bin = _gamma_cat_probs[j];
            inference_process *p_new_process = new inference_process(_ost, _p_lambda, _lambda_multipliers[j], _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _rootdist_vec); // if a single _lambda_multiplier, how do we do it?
            bundle.add(p_new_process);
            
            cout << "  Started inference process " << j+1 << endl;
            
	}

	    _inference_bundles.push_back(bundle);
    }
}

void gamma_core::print_results(std::ostream& ost)
{
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

double base_core::infer_processes()
{
    initialize_rootdist_if_necessary();
    equilibrium_frequency eq(_rootdist_vec);
    std::vector<double> all_families_likelihood(processes.size());
    // prune all the families with the same lambda
    for (int i = 0; i < processes.size(); ++i) {
        auto partial_likelihood = processes[i]->prune();
        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = eq.compute(j);
            full[j] = partial_likelihood[j]*eq_freq;
        	// cout << "Likelihood " << count << ": Partial " << partial_likelihood[j] << ", eq freq: " << eq_freq << ", Full " << full[j] << endl;
        }
        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0);
    }

    double multi = std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 1.0, std::multiplies<double>());

    cout << "Final answer: " << multi << std::endl;

    return multi;
}

//! Infer bundle
double gamma_core::infer_processes() {

    initialize_rootdist_if_necessary();

    equilibrium_frequency eq(_rootdist_vec);
    std::vector<double> all_bundles_likelihood(_inference_bundles.size());

    for (int i = 0; i < _inference_bundles.size(); ++i) {
        cout << endl << "About to prune a gamma bundle." << endl;
        std::vector<family_info_stash> stash = _inference_bundles[i].prune(_gamma_cat_probs, &eq);

        all_bundles_likelihood[i] = 0;
        double family_likelihood = 0;
        for (size_t k = 0; k < stash.size(); ++k)
        {
            stash[k].family_likelihood += stash[k].category_likelihood;
            all_bundles_likelihood[i] = stash[k].family_likelihood;
            stash[k].family_id = i;
        }

        cout << "Likelihood of family " << i << " = " << all_bundles_likelihood[i] << std::endl;
    }

    double multi = std::accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 1.0, std::multiplies<double>());

    cout << "Final answer: " << multi << std::endl;

    return multi;
}

gamma_core::~gamma_core()
{
    for (size_t i = 0; i < _inference_bundles.size(); ++i)
        _inference_bundles[i].clear();
    _inference_bundles.clear();
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

void gamma_bundle::clear()
{
    for (size_t i = 0; i < processes.size(); ++i)
    {
        delete processes[i];
    }
    processes.clear();
}