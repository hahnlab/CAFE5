#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include "core.h"
#include "user_data.h"
#include "matrix_cache.h"
#include "gamma_core.h"
#include "base_model.h"

std::vector<model *> build_models(const input_parameters& user_input, user_data& user_data) {

    model *p_model = NULL;

    std::vector<gene_family> *p_gene_families = &user_data.gene_families;

    if (user_input.is_simulating) {
        p_gene_families = NULL;
    }

    if (user_input.fixed_alpha > 0 || user_input.n_gamma_cats > 1)
    {
        auto gmodel = new gamma_model(user_data.p_lambda, user_data.p_tree, &user_data.gene_families, user_data.max_family_size, user_data.max_root_family_size,
            user_input.n_gamma_cats, user_input.fixed_alpha, user_data.p_error_model);
#ifndef SILENT
        if (user_input.fixed_alpha >= 0)
            gmodel->write_probabilities(cout);
#endif
        p_model = gmodel;
    }
    else
    {
        p_model = new base_model(user_data.p_lambda, user_data.p_tree, p_gene_families, user_data.max_family_size, user_data.max_root_family_size, user_data.p_error_model);
    }

    return std::vector<model *>{p_model};
}

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r)
{
    ost << r.family_id << "\t" << r.lambda_multiplier << "\t" << r.category_likelihood << "\t" << r.family_likelihood;
    ost << "\t" << r.posterior_probability << "\t" << (r.significant ? "*" : "N/S");
    return ost;
}

model::model(lambda* p_lambda,
    const clade *p_tree,
    const vector<gene_family> *p_gene_families,
    int max_family_size,
    int max_root_family_size,
    error_model *p_error_model) :
    _ost(cout), _p_lambda(p_lambda), _p_tree(p_tree), _p_gene_families(p_gene_families), _max_family_size(max_family_size),
    _max_root_family_size(max_root_family_size), _p_error_model(p_error_model) 
{
    if (_p_gene_families)
        references = build_reference_list(*_p_gene_families);
}

std::size_t model::get_gene_family_count() const {
    return _p_gene_families->size();
}

//! Set max family sizes and max root family sizes for INFERENCE
void model::set_max_sizes(int max_family_size, int max_root_family_size) {
    _max_family_size = max_family_size;
    _max_root_family_size = max_root_family_size;
}

class lambda_counter
{
public:
    std::set<int> unique_lambdas;
    void operator()(const clade *p_node)
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
    clademap<double>& _depths;
    double _depth;
public:
    depth_finder(clademap<double>& depths, double depth) : _depths(depths), _depth(depth) {}
    void operator()(const clade *c)
    {
        _depths[c] = _depth + c->get_branch_length();
        depth_finder dp(_depths, _depth + c->get_branch_length());
        c->apply_to_descendants(dp);
    }
};

void model::print_node_depths(std::ostream& ost)
{
    clademap<double> depths;
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

void model::write_vital_statistics(std::ostream& ost, double final_likelihood)
{
    ost << "Model " << name() << " Result: " << final_likelihood << endl;
    ost << "Lambda: " << *get_lambda() << endl;
    if (_p_error_model)
        ost << "Epsilon: " << _p_error_model->get_epsilons()[0] << endl;
}

lambda* model::get_simulation_lambda(const user_data& data)
{
    return data.p_lambda->clone();
}

void branch_length_finder::operator()(const clade *c)
{
    if (c->get_branch_length() > 0.0)
        _result.insert(c->get_branch_length());
}

double branch_length_finder::longest() const
{
    return *max_element(_result.begin(), _result.end());
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by creating a likelihood_computer and calling it on all modes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_prune(const gene_family& gf, matrix_cache& calc, const lambda *_lambda, const clade *_p_tree, double _lambda_multiplier, int _max_root_family_size, int _max_family_size)
{
    unique_ptr<lambda> multiplier(_lambda->multiply(_lambda_multiplier));
    std::map<const clade *, std::vector<double> > _probabilities;
    initialize_probabilities(_p_tree, _probabilities, _max_root_family_size, _max_family_size);
    auto fn = [&](const clade *c) { compute_node_probability(c, gf, NULL, _probabilities, _max_root_family_size, _max_family_size, multiplier.get(), calc); };
    _p_tree->apply_reverse_level_order(fn);

    return _probabilities.at(_p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

void event_monitor::Event_InferenceAttempt_Started() 
{ 
#ifndef SILENT
    std::cout << "-lnL: " << log(0) << std::endl;
#endif
    attempts++;
}

void event_monitor::Event_Reconstruction_Started(std::string model)
{
#ifndef SILENT
    cout << "Starting reconstruction processes for " << model << " model" << endl;
#endif
}

void event_monitor::Event_Reconstruction_Complete()
{
#ifndef SILENT
    cout << "Done!" << endl;
#endif
}

void event_monitor::Event_InferenceAttempt_Complete(double final_likelihood)
{
#ifndef SILENT
    std::cout << "-lnL: " << final_likelihood << std::endl;
#endif
}

void event_monitor::summarize(std::ostream& ost) const
{
    if (attempts == 0)
    {
        ost << "No attempts made\n";
        return;
    }
    ost << this->attempts << " values were attempted (" << round(double(rejects) / double(attempts) * 100) << "% rejected)\n";
    if (!failure_count.empty())
    {
        auto failures = [](const pair<string, int>& a, const pair<string, int>& b) { return a.second < b.second; };
        auto worst_performing_family = std::max_element(failure_count.begin(), failure_count.end(), failures);
        if (worst_performing_family->second * 5 > (attempts - rejects))    // at least one family had 20% rejections
        {
            ost << "The following families had failure rates >20% of the time:\n";
            for (auto& a : this->failure_count)
            {
                if (a.second * 5 > (attempts - rejects))
                    ost << a.first << " had " << a.second << " failures\n";
            }
        }
    }
}

