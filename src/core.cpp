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
#include "error_model.h"

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
        error_model* p_error_model = user_data.p_error_model;
        if (user_input.use_error_model && !p_error_model)
        {
            p_error_model = new error_model();
            p_error_model->set_probabilities(0, { 0, .95, 0.05 });
            p_error_model->set_probabilities(user_data.max_family_size, { 0.05, .9, 0.05 });
        }

        p_model = new base_model(user_data.p_lambda, user_data.p_tree, p_gene_families, user_data.max_family_size, user_data.max_root_family_size, p_error_model);
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

void model::initialize_lambda(clade *p_lambda_tree)
{
    lambda *p_lambda = NULL;
    if (p_lambda_tree != NULL)
    {
        std::set<int> unique_lambdas;
        auto fn = [&unique_lambdas](const clade *p_node) { unique_lambdas.insert(p_node->get_lambda_index()); };
        p_lambda_tree->apply_prefix_order(fn);
        auto node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
        p_lambda = new multiple_lambda(node_name_to_lambda_index, std::vector<double>(unique_lambdas.size()));
        cout << "Searching for " << unique_lambdas.size() << " lambdas" << endl;
    }
    else
    {
        p_lambda = new single_lambda(0.0);
    }

    _p_lambda = p_lambda;
}

void model::write_vital_statistics(std::ostream& ost, double final_likelihood)
{
    ost << "Model " << name() << " Result: " << final_likelihood << endl;
    ost << "Lambda: " << *get_lambda() << endl;
    if (_p_error_model)
        ost << "Epsilon: " << _p_error_model->get_epsilons()[0] << endl;
}

lambda* model::get_simulation_lambda()
{
    return _p_lambda->clone();
}

void model::write_error_model(std::ostream& ost)
{
    auto em = _p_error_model;
    if (!em)
    {
        em = new error_model();
        em->set_probabilities(_max_family_size, { 0, 1, 0 });
    }
    write_error_model_file(ost, *em);
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by calling \ref compute_node_probability on all nodes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_prune(const gene_family& gf, matrix_cache& calc, const lambda *p_lambda, const error_model* p_error_model, const clade *p_tree, double lambda_multiplier, int max_root_family_size, int max_family_size)
{
    unique_ptr<lambda> multiplier(p_lambda->multiply(lambda_multiplier));
    clademap<std::vector<double>> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node].resize(node->is_root() ? max_root_family_size : max_family_size + 1); };
    p_tree->apply_reverse_level_order(init_func);

    auto compute_func = [&](const clade *c) { compute_node_probability(c, gf, p_error_model, probabilities, max_root_family_size, max_family_size, multiplier.get(), calc); };
    p_tree->apply_reverse_level_order(compute_func);

    return probabilities.at(p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

void event_monitor::Event_InferenceAttempt_Started() 
{ 
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
    cout << "Done!\n" << endl;
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

