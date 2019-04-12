#include <cmath>
#include <numeric>
#include <limits>
#include <random>

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

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);

    matrix_cache calc(max(_max_root_family_size, _max_family_size) + 1);
    calc.precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());

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

    branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);

    if (data.p_error_model)
    {
        return new lambda_epsilon_optimizer(this, _p_error_model, data.p_prior.get(), data.rootdist, _p_lambda, finder.longest());
    }
    else
    {
        return new lambda_optimizer(_p_lambda, this, data.p_prior.get(), finder.longest(), data.rootdist);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior)
{
    _monitor.Event_Reconstruction_Started("Base");

    std::vector<clademap<int>> reconstructed_states;
    std::vector<clademap<family_size_change>> increase_decrease_map;
    std::vector<string> family_ids;

    family_ids.resize(_p_gene_families->size());
    transform(_p_gene_families->begin(), _p_gene_families->end(), family_ids.begin(), [](const gene_family& fam) { return fam.id(); });

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    p_calc->precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());

    reconstructed_states.resize(_p_gene_families->size());
    increase_decrease_map.resize(_p_gene_families->size());
    for (size_t i = 0; i<_p_gene_families->size(); ++i)
    {
        reconstruct_gene_families(_p_lambda, _p_tree, _max_family_size, _max_root_family_size,
            &_p_gene_families->at(i), p_calc, p_prior, reconstructed_states[i]);
        compute_increase_decrease(reconstructed_states[i], increase_decrease_map[i]);
    }

    _monitor.Event_Reconstruction_Complete();
    return new base_model_reconstruction(reconstructed_states, increase_decrease_map, family_ids);
}

void base_model::prepare_matrices_for_simulation(matrix_cache& cache)
{
    unique_ptr<lambda> perturbed_lambda(_p_lambda->multiply(simulation_lambda_multiplier));
    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    cache.precalculate_matrices(get_lambda_values(perturbed_lambda.get()), lengths.result());
}

lambda* base_model::get_simulation_lambda(const user_data& data)
{
    return data.p_lambda->multiply(simulation_lambda_multiplier);
}

void base_model::perturb_lambda()
{
    normal_distribution<double> dist(1.0, 0.3);
    simulation_lambda_multiplier = dist(randomizer_engine);
}

void base_model_reconstruction::print_reconstructed_states(std::ostream& ost, const std::vector<gene_family>& gene_families, const clade *p_tree) {
    if (_reconstructed_family_counts.empty())
        return;

    auto order = p_tree->find_internal_nodes();

    ost << "#NEXUS\nBEGIN TREES;\n";
    for (size_t i = 0; i<gene_families.size(); ++i)
    {
        auto& gene_family = gene_families[i];
        auto g = [i, gene_family, this](const clade *node) {
            int value = node->is_leaf() ? gene_family.get_species_size(node->get_taxon_name()) : _reconstructed_family_counts[i].at(node);
            return to_string(value);
        };

        auto f = [g, order, this](const clade *node) {
            return newick_node(node, order, g);
        };

        ost << "  TREE " << gene_family.id() << " = ";
        p_tree->write_newick(ost, f);

        ost << ';' << endl;
    }
    ost << "END;\n";
}

increase_decrease get_increases_decreases(const clademap<family_size_change>& increase_decrease_map, const cladevector& order, double pvalue, string family_id)
{
    increase_decrease result;
    result.change.resize(order.size());
    result.gene_family_id = family_id;
    result.pvalue = pvalue;

    transform(order.begin(), order.end(), result.change.begin(), [increase_decrease_map](const clade *taxon)->family_size_change {
        return increase_decrease_map.at(taxon);
    });

    return result;
}

void base_model_reconstruction::print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, const std::vector<double>& pvalues) {
    if (_reconstructed_family_counts.size() != pvalues.size())
    {
        throw std::runtime_error("No pvalues found for family");
    }
    if (_family_increase_decrease.empty())
    {
        ost << "No increases or decreases recorded\n";
        return;
    }

    ost << "#FamilyID\tpvalue\t*\t";
    for (auto& it : order) {
        ost << it->get_taxon_name() << "\t";
    }
    ost << endl;

    for (size_t i = 0; i < _family_increase_decrease.size(); ++i) {
        ost << get_increases_decreases(_family_increase_decrease[i], order, pvalues[i], _family_ids[i]);
    }
}

void base_model_reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order) {
    if (_reconstructed_family_counts.empty())
    {
        ost << "No increases or decreases recorded\n";
        return;
    }

    clademap<pair<int, int>> items;

    for (size_t i = 0; i < _family_increase_decrease.size(); ++i) {
        auto incdec = get_increases_decreases(_family_increase_decrease[i], order, 0.0, _family_ids[i]);
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (incdec.change[i] == Increase)
                items[order[i]].first++;
            if (incdec.change[i] == Decrease)
                items[order[i]].second++;
        }
    }

    ost << "#Taxon_ID\tIncrease/Decrease\n";
    for (auto& it : items) {
        ost << it.first->get_taxon_name() << "\t";
        ost << it.second.first << "/" << it.second.second << endl;
    }
}

