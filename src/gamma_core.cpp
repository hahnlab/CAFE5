#include <assert.h>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <random>
#include <sstream>

#include "gamma_core.h"
#include "gamma.h"
#include "root_equilibrium_distribution.h"
#include "gene_family_reconstructor.h"
#include "matrix_cache.h"
#include "gene_family.h"
#include "user_data.h"
#include "optimizer_scorer.h"
#include "root_distribution.h"
#include "simulator.h"

extern mt19937 randomizer_engine;

gamma_model::gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
    int max_root_family_size, int n_gamma_cats, double fixed_alpha, error_model* p_error_model) :
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_error_model) {

    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
    if (p_gene_families)
        _category_likelihoods.resize(p_gene_families->size());
    set_alpha(fixed_alpha);
}

gamma_model::gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
    int max_root_family_size, std::vector<double> gamma_categories, std::vector<double> multipliers, error_model *p_error_model) :
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_error_model)
{
    _gamma_cat_probs = gamma_categories;
    _lambda_multipliers = multipliers;
}

void gamma_model::write_vital_statistics(std::ostream& ost, double final_likelihood)
{
    model::write_vital_statistics(ost, final_likelihood);
    ost << "Alpha: " << _alpha << endl;
}

void gamma_model::write_family_likelihoods(std::ostream& ost)
{
    ost << "#FamilyID\tGamma Cat Median\tLikelihood of Category\tLikelihood of Family\tPosterior Probability\tSignificant" << endl;

    std::ostream_iterator<family_info_stash> out_it(ost, "\n");
    std::copy(results.begin(), results.end(), out_it);
}

//! Set alpha for gamma distribution
void gamma_model::set_alpha(double alpha) {

    _alpha = alpha;
    if (_gamma_cat_probs.size() > 1)
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference

}

string comma_separated(const std::vector<double>& items)
{
    string s;
    for (auto i : items)
        s += (s.empty() ? "" : ",") + to_string(i);
    return s;
}

void gamma_model::write_probabilities(ostream& ost)
{
    ost << "Gamma cat probs are: " << comma_separated(_gamma_cat_probs) << endl;
    ost << "Lambda multipliers are: " << comma_separated(_lambda_multipliers) << endl;
}

lambda* gamma_model::get_simulation_lambda(const user_data& data)
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());

    return data.p_lambda->multiply(_lambda_multipliers[dist(randomizer_engine)]);
}

std::vector<double> gamma_model::get_posterior_probabilities(std::vector<double> cat_likelihoods)
{
    size_t process_count = cat_likelihoods.size();

    vector<double> numerators(process_count);
    transform(cat_likelihoods.begin(), cat_likelihoods.end(), _gamma_cat_probs.begin(), numerators.begin(), multiplies<double>());

    double denominator = accumulate(numerators.begin(), numerators.end(), 0.0);
    vector<double> posterior_probabilities(process_count);
    transform(numerators.begin(), numerators.end(), posterior_probabilities.begin(), bind2nd(divides<double>(), denominator));

    return posterior_probabilities;
}

void gamma_model::prepare_matrices_for_simulation(matrix_cache& cache)
{
    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    vector<double> multipliers;
    for (auto multiplier : _lambda_multipliers)
    {
        unique_ptr<lambda> mult(_p_lambda->multiply(multiplier));
        auto values = get_lambda_values(mult.get());
        multipliers.insert(multipliers.end(), values.begin(), values.end());
    }
    cache.precalculate_matrices(multipliers, lengths.result());
}

void gamma_model::perturb_lambda()
{
    if (_gamma_cat_probs.size() == 1)
    {
        // no user cluster value was specified. Select a multiplier at random from the gamma distribution with the given alpha
        gamma_distribution<double> dist(_alpha, 1 / _alpha);
        _lambda_multipliers[0] = dist(randomizer_engine);
        _gamma_cat_probs[0] = 1;
    }
    else
    {
        // select multipliers based on the clusters, modifying the actual lambda selected
        // by a normal distribution based around the multipliers selected from the gamma
        // distribution with the given alpha
        
        // first, reset the multipliers back to their initial values based on the alpha
        get_gamma(_gamma_cat_probs, _lambda_multipliers, _alpha);

        auto new_multipliers = _lambda_multipliers;
        for (size_t i = 0; i < _lambda_multipliers.size(); ++i)
        {
            double stddev;
            if (i == 0)
            {
                stddev = _lambda_multipliers[0] / 3.0;
            }
            else if (i == _lambda_multipliers.size() - 1)
            {
                stddev = (_lambda_multipliers[i] - _lambda_multipliers[i - 1]) / 3.0;
            }
            else
            {
                stddev = (_lambda_multipliers[i + 1] - _lambda_multipliers[i - 1]) / 6.0;
            }
            normal_distribution<double> dist(_lambda_multipliers[i], stddev);
            new_multipliers[i] = dist(randomizer_engine);
        }
        _lambda_multipliers.swap(new_multipliers);
    }

#ifndef SILENT
    write_probabilities(cout);
#endif
}

bool gamma_model::can_infer() const
{
    if (!_p_lambda->is_valid())
        return false;

    if (_alpha < 0)
        return false;

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    auto v = get_lambda_values(_p_lambda);

    double longest_branch = lengths.longest();
    double largest_multiplier = *max_element(_lambda_multipliers.begin(), _lambda_multipliers.end());
    double largest_lambda = *max_element(v.begin(), v.end());

    if (matrix_cache::is_saturated(longest_branch, largest_multiplier*largest_lambda))
        return false;

    return true;
}

bool gamma_model::prune(const gene_family& family, root_equilibrium_distribution *eq, matrix_cache& calc, const lambda *p_lambda,
    std::vector<double>& category_likelihoods) 
{
    category_likelihoods.clear();

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        auto partial_likelihood = inference_prune(family, calc, p_lambda, _p_tree, _lambda_multipliers[k], _max_root_family_size, _max_family_size);
        if (accumulate(partial_likelihood.begin(), partial_likelihood.end(), 0.0) == 0.0)
            return false;   // saturation

        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = eq->compute(j);
            full[j] = partial_likelihood[j] * eq_freq;
        }

        //        _category_likelihoods.push_back(accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k]); // sum over all sizes (Felsenstein's approach)
        category_likelihoods.push_back(*max_element(full.begin(), full.end()) * _gamma_cat_probs[k]); // get max (CAFE's approach)
    }

    return true;
}

//! Infer bundle
double gamma_model::infer_family_likelihoods(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map, const lambda *p_lambda) {

    _monitor.Event_InferenceAttempt_Started();

    if (!can_infer())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    using namespace std;
    root_distribution rd;
    if (root_distribution_map.size() > 0)
    {
        rd.vectorize(root_distribution_map);
    }
    else
    {
        rd.vectorize_uniform(_max_root_family_size);
    }

    prior->initialize(&rd);
    vector<double> all_bundles_likelihood(_p_gene_families->size());

    vector<bool> failure(_p_gene_families->size());
    matrix_cache calc(max(_max_root_family_size, _max_family_size) + 1);
    prepare_matrices_for_simulation(calc);

    vector<vector<family_info_stash>> pruning_results(_p_gene_families->size());

#pragma omp parallel for
    for (size_t i = 0; i < _p_gene_families->size(); i++) {
        auto& cat_likelihoods = _category_likelihoods[i];

        if (prune(_p_gene_families->at(i), prior, calc, p_lambda, cat_likelihoods))
        {
            double family_likelihood = accumulate(cat_likelihoods.begin(), cat_likelihoods.end(), 0.0);

            vector<double> posterior_probabilities = get_posterior_probabilities(cat_likelihoods);

            pruning_results[i].resize(cat_likelihoods.size());
            for (size_t k = 0; k < cat_likelihoods.size(); ++k)
            {
                pruning_results[i][k] = family_info_stash(_p_gene_families->at(i).id(),_lambda_multipliers[k], cat_likelihoods[k],
                    family_likelihood, posterior_probabilities[k], posterior_probabilities[k] > 0.95);
                //            cout << "Bundle " << i << " Process " << k << " family likelihood = " << family_likelihood << endl;
            }
            all_bundles_likelihood[i] = std::log(family_likelihood);
        }
        else
        {
            // we got here because one of the gamma categories was saturated - reject this 
            failure[i] = true;
        }
    }

    if (find(failure.begin(), failure.end(), true) != failure.end())
    {
        for (size_t i = 0; i < _p_gene_families->size(); i++) {
            if (failure[i])
            {
                _monitor.Event_InferenceAttempt_Saturation(_p_gene_families->at(i).id());
            }
        }
        return -log(0);
    }
    for (auto& stashes : pruning_results)
    {
        for (auto& stash : stashes)
        {
            results.push_back(stash);
        }
    }
    double final_likelihood = -accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 0.0);

    _monitor.Event_InferenceAttempt_Complete(final_likelihood);
    return final_likelihood;
}

inference_optimizer_scorer *gamma_model::get_lambda_optimizer(user_data& data)
{
    bool estimate_lambda = data.p_lambda == NULL;
    bool estimate_alpha = _alpha <= 0.0;

    if (estimate_lambda && estimate_alpha)
    {
        branch_length_finder finder;
        _p_tree->apply_prefix_order(finder);

        initialize_lambda(data.p_lambda_tree);
        return new gamma_lambda_optimizer(_p_lambda, this, data.p_prior.get(), data.rootdist, finder.longest());
    }
    else if (estimate_lambda && !estimate_alpha)
    {
        branch_length_finder finder;
        _p_tree->apply_prefix_order(finder);

        initialize_lambda(data.p_lambda_tree);
        return new lambda_optimizer(_p_lambda, this, data.p_prior.get(), finder.longest(), data.rootdist);
    }
    else if (!estimate_lambda && estimate_alpha)
    {
        _p_lambda = data.p_lambda->clone();
        return new gamma_optimizer(this, data.p_prior.get(), data.rootdist);
    }
    else
    {
        return nullptr;
    }
}

clademap<double> get_weighted_averages(const std::vector<reconstructed_family<int>>& m, const vector<double>& probabilities)
{
    cladevector nodes(m[0].clade_counts.size());
    std::transform(m[0].clade_counts.begin(), m[0].clade_counts.end(), nodes.begin(), [](std::pair<const clade *, int> v) { return v.first;  });

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<probabilities.size(); ++i)
        {
            val += probabilities[i] * double(m[i].clade_counts.at(node));
        }
        result[node] = val;
    }

    return result;
}

void gamma_model::reconstruct_family(const gene_family& family, matrix_cache *calc, root_equilibrium_distribution*prior, gamma_model_reconstruction::gamma_reconstruction& rc) const
{
    auto& cat_rec = rc.category_reconstruction;

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        unique_ptr<lambda> ml(_p_lambda->multiply(_lambda_multipliers[k]));
        reconstruct_gene_families(ml.get(), _p_tree, _max_family_size, _max_root_family_size, &family, calc, prior, cat_rec[k].clade_counts);
        compute_increase_decrease(cat_rec[k].clade_counts, cat_rec[k].size_deltas);
    }

    // multiply every reconstruction by gamma_cat_prob
    rc.reconstruction.clade_counts = get_weighted_averages(cat_rec, _gamma_cat_probs);

    compute_increase_decrease(rc.reconstruction.clade_counts, rc.reconstruction.size_deltas);
}

reconstruction* gamma_model::reconstruct_ancestral_states(matrix_cache *calc, root_equilibrium_distribution*prior)
{
    _monitor.Event_Reconstruction_Started("Gamma");
    gamma_model_reconstruction* result = new gamma_model_reconstruction(_p_gene_families->size(), _lambda_multipliers);

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    auto values = get_lambda_values(_p_lambda);
    vector<double> all;
    for (double multiplier : _lambda_multipliers)
    {
        for (double lambda : values)
        {
            all.push_back(lambda*multiplier);
        }
    }

    calc->precalculate_matrices(all, lengths.result());

#pragma omp parallel for
    for (size_t i = 0; i<_p_gene_families->size(); ++i)
    {
        result->_families[i].reconstruction.id = _p_gene_families->at(i).id();

        reconstruct_family(_p_gene_families->at(i), calc, prior, result->_families[i]);

        result->_families[i]._category_likelihoods = _category_likelihoods[i];
    }

    _monitor.Event_Reconstruction_Complete();

    return result;
}

void gamma_model_reconstruction::print_reconstructed_states(std::ostream& ost, const cladevector& order, const std::vector<gene_family>& gene_families, const clade *p_tree)
{
    if (_families.empty())
        return;

    auto rec = _families[0];

    ost << "#NEXUS\nBEGIN TREES;\n";
    for (size_t i = 0; i<gene_families.size(); ++i)
    {
        auto& gene_family = gene_families[i];

        auto g = [i, gene_family, this](const clade *node) {
            std::ostringstream ost;

            if (node->is_leaf())
            {
                ost << gene_family.get_species_size(node->get_taxon_name());
            }
            else
            {
                for (auto& r : _families[i].category_reconstruction)
                {
                    ost << r.clade_counts.at(node) << '_';
                }
                ost << std::round(_families[i].reconstruction.clade_counts.at(node));
            }
            return ost.str();
        };

        auto f = [order, g, this](const clade *node) {
            return newick_node(node, order, g);
        };

        ost << "  TREE " << gene_family.id() << " = ";
        p_tree->write_newick(ost, f);

        ost << ';' << endl;
    }
    ost << "END;\n\n";

    ost << "BEGIN LAMBDA_MULTIPLIERS;\n";
    for (auto& lm : _lambda_multipliers)
    {
        ost << "  " << lm << ";\n";
    }
    ost << "END;\n";
    ost << endl;
}

increase_decrease get_increases_decreases(const gamma_model_reconstruction::gamma_reconstruction& rc, const cladevector& order, double pvalue)
{
    increase_decrease result;
    result.change.resize(order.size());
    result.gene_family_id = rc.reconstruction.id;
    result.pvalue = pvalue;

    transform(order.begin(), order.end(), result.change.begin(), [rc](const clade *taxon)->family_size_change {
        if (taxon->is_leaf() || taxon->is_root())
            return Constant;
        else
            return rc.reconstruction.size_deltas.at(taxon);
    });

    result.category_likelihoods = rc._category_likelihoods;
    return result;
}


void gamma_model_reconstruction::print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, const std::vector<double>& pvalues)
{
    reconstruction::print_increases_decreases_by_family(ost, order, pvalues, _families.size(), [this, order, pvalues](int family_index) {
        return get_increases_decreases(_families[family_index], order, pvalues[family_index]);
        });
}

void gamma_model_reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order)
{
    reconstruction::print_increases_decreases_by_clade(ost, order, _families.size(), [this, order](int family_index) {
        return get_increases_decreases(_families[family_index], order, 0.0);
        });
}

