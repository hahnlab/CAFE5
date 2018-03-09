#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include "clade.h"
#include "core.h"
#include "base_model.h"
#include "gamma_core.h"
#include "family_generator.h"
#include "process.h"
#include "poisson.h"
#include "root_equilibrium_distribution.h"
#include "reconstruction_process.h"

std::vector<model *> build_models(const input_parameters& my_input_parameters, 
    clade *p_tree, 
    lambda *p_lambda, 
    std::vector<gene_family>* p_gene_families, 
    int max_family_size, 
    int max_root_family_size,
    error_model *p_error_model) {
    std::vector<model *> models;
    
    /* If estimating or computing (-i; or -i + -l) and not simulating */
    if (my_input_parameters.nsims == 0 && !my_input_parameters.input_file_path.empty()) {
        
        /* Base core is always used (in both computation and estimation) */
        models.push_back(new base_model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, NULL, p_error_model));
        
        /* Gamma core is only used in estimation */
        if (p_lambda == NULL && my_input_parameters.n_gamma_cats > 1) {
            auto model = new gamma_model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size,
                my_input_parameters.n_gamma_cats, my_input_parameters.fixed_alpha, NULL, p_error_model);
            model->write_probabilities(cout);
            models.push_back(model);
        }
    }
    
    /* If simulating (-s) */
    else if (my_input_parameters.nsims > 0 || !my_input_parameters.rootdist.empty()) {
        
        /* -f */
        std::map<int, int> *p_rootdist_map = read_rootdist(my_input_parameters.rootdist); // in map form

        // Either use base core or gamma core when simulating
        if (my_input_parameters.n_gamma_cats > 1) { models.push_back(new gamma_model(p_lambda, p_tree, NULL, 
            max_family_size, max_root_family_size, my_input_parameters.n_gamma_cats, my_input_parameters.fixed_alpha,
            p_rootdist_map, p_error_model)); }
        else { models.push_back(new base_model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_rootdist_map, p_error_model)); }
    }

    return models;
}

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r)
{
    ost << r.family_id << "\t" << r.lambda_multiplier << "\t" << r.category_likelihood << "\t" << r.family_likelihood;
    ost << "\t" << r.posterior_probability << "\t" << (r.significant ? "*" : "N/S");
    return ost;
}

model::model(lambda* p_lambda,
    clade *p_tree,
    vector<gene_family> *p_gene_families,
    int max_family_size,
    int max_root_family_size,
    error_model *p_error_model) :
    _ost(cout), _p_lambda(p_lambda), _p_tree(p_tree), _p_gene_families(p_gene_families), _max_family_size(max_family_size),
    _max_root_family_size(max_root_family_size), _p_error_model(p_error_model) 
{
    if (_p_gene_families)
        references = build_reference_list(*_p_gene_families);
}

//! Set pointer to lambda in core class
void model::set_tree(clade *p_tree) {
    _p_tree = p_tree;
}

//! Set pointer to vector of gene family class instances
void model::set_gene_families(std::vector<gene_family> *p_gene_families) {
    _p_gene_families = p_gene_families;
}

//! Set max family sizes and max root family sizes for INFERENCE
void model::set_max_sizes(int max_family_size, int max_root_family_size) {
    _max_family_size = max_family_size;
    _max_root_family_size = max_root_family_size;
}

//! Set root distribution vector
void model::set_rootdist_vec(std::vector<int> rootdist_vec) {
    _rootdist_vec = rootdist_vec;
}

//! Set total number of families to simulate
void model::set_total_n_families_sim(int total_n_families_sim) {
    _total_n_families_sim = total_n_families_sim;
}

void model::start_sim_processes() {

    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes.push_back(create_simulation_process(i));
    }

    // cout << _sim_processes.size() << " processes have been started." << endl;
}

//! Run simulations in all processes, in series... (TODO: in parallel!)
void model::simulate_processes() {
    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->run_simulation();
    }
}

void model::print_processes(std::ostream& ost) {

    if (_sim_processes.empty())
    {
        cerr << "No simulations created" << endl;
        return;
    }
    trial *sim = _sim_processes[0]->get_simulation();
    for (trial::iterator it = sim->begin(); it != sim->end(); ++it) {
          ost << "#" << it->first->get_taxon_name() << "\n";
    }

    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->print_simulation(ost, i);
    }
}

void model::initialize_rootdist_if_necessary()
{
    if (_rootdist_vec.empty())
    {
        _rootdist_vec.resize(_max_root_family_size);
        std::fill(_rootdist_vec.begin(), _rootdist_vec.end(), 1);
    }

}

//! Print processes' simulations
void model::adjust_family(std::ostream& ost) {
    
    // Printing header
    for (trial::iterator it = _sim_processes[0]->get_simulation()->begin(); it != _sim_processes[0]->get_simulation()->end(); ++it) {
	ost << "#" << it->first->get_taxon_name() << endl;
    }
    
    for (int i = 0; i < _total_n_families_sim; ++i) {
        _sim_processes[i]->print_simulation(cout, i);
    }
}

/* TODO: later this will become a member of the core class, which is the wrapper of the process class*/
void model::print_parameter_values() {
    
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

class lambda_counter
{
public:
    std::set<int> unique_lambdas;
    void operator()(clade *p_node)
    {
        unique_lambdas.insert(p_node->get_lambda_index());
    }
};


void model::initialize_lambda(clade *p_lambda_tree)
{
    lambda *p_lambda = NULL;
    if (p_lambda_tree != NULL)
    {
        lambda_counter counter;
        p_lambda_tree->apply_prefix_order(counter);
        auto node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
        p_lambda = new multiple_lambda(node_name_to_lambda_index, std::vector<double>(counter.unique_lambdas.size()));
        cout << "Searching for " << counter.unique_lambdas.size() << " lambdas" << endl;
    }
    else
    {
        p_lambda = new single_lambda(0.0);
    }

    _p_lambda = p_lambda;
}

class depth_finder
{
    std::map<clade *, double>& _depths;
    double _depth;
public:
    depth_finder(std::map<clade *, double>& depths, double depth) : _depths(depths), _depth(depth) {}
    void operator()(clade *c)
    {
        _depths[c] = _depth + c->get_branch_length();
        depth_finder dp(_depths, _depth + c->get_branch_length());
        c->apply_to_descendants(dp);
    }
};

void model::print_node_depths(std::ostream& ost)
{
    std::map<clade *, double> depths;
    depth_finder dp(depths, 0);
    dp(_p_tree);
    
    double max_depth = 0;
    for (auto& x : depths)
    {
        if (x.second > max_depth)
            max_depth = x.second;
    }

    for (auto& x : depths)
    {
        if (!x.first->is_leaf())
        {
            ost << x.first->get_taxon_name() << "\t" << (max_depth - x.second) << endl;
        }
    }
}


void branch_length_finder::operator()(clade *c)
{
    _result.insert(c->get_branch_length());
}

double branch_length_finder::longest() const
{
    return *max_element(_result.begin(), _result.end());
}
