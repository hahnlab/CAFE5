#include <cmath>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>

#include "base_model.h"
#include "gene_family_reconstructor.h"
#include "matrix_cache.h"
#include "gene_family.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "optimizer_scorer.h"
#include "root_distribution.h"
#include "simulator.h"

extern mt19937 randomizer_engine;

base_model::base_model(lambda* p_lambda, const clade *p_tree, const vector<gene_family>* p_gene_families,
    int max_family_size, int max_root_family_size, error_model *p_error_model) :
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_error_model)
{

}

vector<size_t> build_reference_list(const vector<gene_family>& families)
{
    const size_t invalid_index = std::numeric_limits<size_t>::max();

    vector<size_t> reff;
    const int num_families = families.size();
    reff.resize(num_families, invalid_index);
    for (int i = 0; i < num_families; ++i) {
        if (reff[i] != invalid_index) continue;

        reff[i] = i;

        for (int j = i + 1; j < num_families; ++j) {
            if (reff[j] == invalid_index)
            {
                if (families[i].species_size_match(families[j]))
                {
                    reff[j] = i;
                }
            }
        }
    }

    return reff;
}

double base_model::infer_family_likelihoods(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map, const lambda *p_lambda) {
    _monitor.Event_InferenceAttempt_Started();

    if (!_p_lambda->is_valid())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    root_distribution rd;
    if (root_distribution_map.size() > 0)
    {
        rd.vectorize(root_distribution_map);
    }
    else
    {
        rd.vectorize_uniform(_max_root_family_size);
    }
//    initialize_rootdist_if_necessary();
    prior->initialize(&rd);

    results.resize(_p_gene_families->size());
    std::vector<double> all_families_likelihood(_p_gene_families->size());

    matrix_cache calc(max(_max_root_family_size, _max_family_size) + 1);
    calc.precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());

    vector<vector<double>> partial_likelihoods(_p_gene_families->size());
#pragma omp parallel for
    for (size_t i = 0; i < _p_gene_families->size(); ++i) {
        if (references[i] == i)
            partial_likelihoods[i] = inference_prune(_p_gene_families->at(i), calc, _p_lambda, _p_tree, 1.0, _max_root_family_size, _max_family_size);
            // probabilities of various family sizes
    }

    // prune all the families with the same lambda
#pragma omp parallel for
    for (size_t i = 0; i < _p_gene_families->size(); ++i) {

        auto& partial_likelihood = partial_likelihoods[references[i]];
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = prior->compute(j);

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }

        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
                                                                             // cout << i << " contribution " << scientific << all_families_likelihood[i] << endl;
        
        results[i] = family_info_stash(_p_gene_families->at(i).id(), 0.0, 0.0, 0.0, all_families_likelihood[i], false);
    }
    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    _monitor.Event_InferenceAttempt_Complete(final_likelihood);

    return final_likelihood;
}

void base_model::write_family_likelihoods(std::ostream& ost)
{
    ost << "#FamilyID\tLikelihood of Family" << endl;
    for (const auto& r : results)
    {
        ost << r.family_id << "\t" << r.posterior_probability << endl;
    }
}

inference_optimizer_scorer *base_model::get_lambda_optimizer(user_data& data)
{
    if (data.p_lambda != NULL)  // already have a lambda, nothing we want to optimize
        return nullptr;

    initialize_lambda(data.p_lambda_tree);

    auto lengths = _p_tree->get_branch_lengths();
    auto longest_branch = *max_element(lengths.begin(), lengths.end());

    if (data.p_error_model)
    {
        return new lambda_epsilon_optimizer(this, _p_error_model, data.p_prior.get(), data.rootdist, _p_lambda, longest_branch);
    }
    else
    {
        return new lambda_optimizer(_p_lambda, this, data.p_prior.get(), longest_branch, data.rootdist);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(const vector<const gene_family*>& families, matrix_cache *p_calc, root_equilibrium_distribution* p_prior)
{
    _monitor.Event_Reconstruction_Started("Base");

    auto result = new base_model_reconstruction();

    p_calc->precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());

    for (size_t i = 0; i< families.size(); ++i)
    {
        reconstruct_gene_family(_p_lambda, _p_tree, _max_family_size, _max_root_family_size,
            families[i], p_calc, p_prior, result->_reconstructions[families[i]->id()]);
    }

    _monitor.Event_Reconstruction_Complete();

    return result;
}

void base_model::prepare_matrices_for_simulation(matrix_cache& cache)
{
    unique_ptr<lambda> perturbed_lambda(get_simulation_lambda());
    cache.precalculate_matrices(get_lambda_values(perturbed_lambda.get()), _p_tree->get_branch_lengths());
}

lambda* base_model::get_simulation_lambda()
{
    return _p_lambda->multiply(simulation_lambda_multiplier);
}

void base_model::perturb_lambda()
{
    normal_distribution<double> dist(1.0, 0.3);
    simulation_lambda_multiplier = dist(randomizer_engine);
}

std::string base_model_reconstruction::get_reconstructed_state(const gene_family&gf, const clade* node)
{
    int value = node->is_leaf() ? gf.get_species_size(node->get_taxon_name()) : _reconstructions[gf.id()].at(node);
    return to_string(value);
}

int base_model_reconstruction::get_difference_from_parent(const gene_family* gf, const clade* c)
{
    if (c->is_root())
        return 0;
    int val = c->is_leaf() ? gf->get_species_size(c->get_taxon_name()) : _reconstructions[gf->id()].at(c);
    int parent_val = _reconstructions[gf->id()].at(c->get_parent());

    return val - parent_val;
}

int base_model_reconstruction::get_node_count(const gene_family& gf, const clade *c)
{
    return _reconstructions[gf.id()].at(c);
}

int base_model_reconstruction::reconstructed_size(const gene_family& family, const clade* clade) const
{
    if (clade->is_leaf())
        return family.get_species_size(clade->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw std::runtime_error("Family " + family.id() + " was not reconstructed");
    auto& c = _reconstructions.at(family.id());
    if (c.find(clade) == c.end())
        throw std::runtime_error("Clade '" + clade->get_taxon_name() + "' was not reconstructed for family " + family.id());

    return _reconstructions.at(family.id()).at(clade);
}
