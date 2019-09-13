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

lambda* gamma_model::get_simulation_lambda()
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());

    return _p_lambda->multiply(_lambda_multipliers[dist(randomizer_engine)]);
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
    vector<double> multipliers;
    for (auto multiplier : _lambda_multipliers)
    {
        unique_ptr<lambda> mult(_p_lambda->multiply(multiplier));
        auto values = get_lambda_values(mult.get());
        multipliers.insert(multipliers.end(), values.begin(), values.end());
    }
    cache.precalculate_matrices(multipliers, _p_tree->get_branch_lengths());
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

    auto v = get_lambda_values(_p_lambda);

    auto lengths = _p_tree->get_branch_lengths();
    auto longest_branch = *max_element(lengths.begin(), lengths.end());
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

    results.clear();

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
        auto lengths = _p_tree->get_branch_lengths();
        auto longest_branch = *max_element(lengths.begin(), lengths.end());

        initialize_lambda(data.p_lambda_tree);
        return new gamma_lambda_optimizer(_p_lambda, this, data.p_prior.get(), data.rootdist, longest_branch);
    }
    else if (estimate_lambda && !estimate_alpha)
    {
        auto lengths = _p_tree->get_branch_lengths();
        auto longest_branch = *max_element(lengths.begin(), lengths.end());

        initialize_lambda(data.p_lambda_tree);
        return new lambda_optimizer(_p_lambda, this, data.p_prior.get(), longest_branch, data.rootdist);
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

clademap<double> get_weighted_averages(const std::vector<clademap<int>>& m, const vector<double>& probabilities)
{
    cladevector nodes(m[0].size());
    std::transform(m[0].begin(), m[0].end(), nodes.begin(), [](std::pair<const clade *, int> v) { return v.first;  });

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<probabilities.size(); ++i)
        {
            val += probabilities[i] * double(m[i].at(node));
        }
        result[node] = val;
    }

    return result;
}

reconstruction* gamma_model::reconstruct_ancestral_states(const vector<const gene_family*>& families, matrix_cache *calc, root_equilibrium_distribution*prior)
{
    _monitor.Event_Reconstruction_Started("Gamma");

    auto values = get_lambda_values(_p_lambda);
    vector<double> all;
    for (double multiplier : _lambda_multipliers)
    {
        for (double lambda : values)
        {
            all.push_back(lambda*multiplier);
        }
    }

    calc->precalculate_matrices(all, _p_tree->get_branch_lengths());

    gamma_model_reconstruction* result = new gamma_model_reconstruction(_lambda_multipliers);
    vector<gamma_model_reconstruction::gamma_reconstruction *> recs(families.size());
    for (size_t i = 0; i < families.size(); ++i)
    {
        recs[i] = &result->_reconstructions[families[i]->id()];
        result->_reconstructions[families[i]->id()]._category_likelihoods = _category_likelihoods[i];
        result->_reconstructions[families[i]->id()].category_reconstruction.resize(_lambda_multipliers.size());
    }


    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        unique_ptr<lambda> ml(_p_lambda->multiply(_lambda_multipliers[k]));

#pragma omp parallel for
        for (size_t i = 0; i < families.size(); ++i)
        {
            reconstruct_gene_family(ml.get(), _p_tree, _max_family_size, _max_root_family_size, families[i], calc, prior, recs[i]->category_reconstruction[k]);
        }
    }

    for (auto reconstruction : recs)
    {
        // multiply every reconstruction by gamma_cat_prob
        reconstruction->reconstruction = get_weighted_averages(reconstruction->category_reconstruction, _gamma_cat_probs);
    }

    _monitor.Event_Reconstruction_Complete();

    return result;
}

bool gamma_model::should_calculate_pvalue(const gene_family& gf) const
{
    vector<const family_info_stash*> values;
    for (auto& s : results) {
        if (s.family_id == gf.id())
            values.push_back(&s);
    }
    auto highest_likelihood = max_element(values.begin(), values.end(), [](const family_info_stash * a, const family_info_stash* b) { return a->category_likelihood < b->category_likelihood; });
    auto fastest_rate = max_element(values.begin(), values.end(), [](const family_info_stash* a, const family_info_stash* b) { return a->lambda_multiplier < b->lambda_multiplier; });
    return fastest_rate == highest_likelihood;
}

lambda* gamma_model::get_pvalue_lambda() const
{
    return _p_lambda->multiply(_lambda_multipliers.back());
}

std::string gamma_model_reconstruction::get_reconstructed_state(const gene_family& gf, const clade* node)
{
    std::ostringstream ost;

    if (node->is_leaf())
    {
        ost << gf.get_species_size(node->get_taxon_name());
    }
    else
    {
        for (auto& r : _reconstructions[gf.id()].category_reconstruction)
        {
            ost << r.at(node) << '_';
        }
        ost << std::round(_reconstructions[gf.id()].reconstruction.at(node));
    }
    return ost.str();
}

void gamma_model_reconstruction::write_nexus_extensions(std::ostream& ost)
{
    ost << "BEGIN LAMBDA_MULTIPLIERS;\n";
    for (auto& lm : _lambda_multipliers)
    {
        ost << "  " << lm << ";\n";
    }
    ost << "END;\n\n";
}

char gamma_model_reconstruction::get_increase_decrease(const gene_family* gf, const clade* c)
{
    clademap<int> size_deltas;
    compute_increase_decrease(_reconstructions[gf->id()].reconstruction, size_deltas);

    int val = 0;
    if (!c->is_leaf() && !c->is_root())
        val = size_deltas.at(c);
    if (val < 0)
        return 'd';
    else if (val > 0)
        return 'i';
    else
        return 'c';
}

int gamma_model_reconstruction::get_delta(const gene_family* gf, const clade* c)
{
    clademap<int> size_deltas;
    compute_increase_decrease(_reconstructions[gf->id()].reconstruction, size_deltas);
    auto& d = size_deltas;
    auto it = d.find(c);
    if (it == d.end())
        return 0;
    return it->second;
}

int gamma_model_reconstruction::get_node_count(const gene_family& gf, const clade* c)
{
    return int(std::round(_reconstructions[gf.id()].reconstruction.at(c)));
}

clademap<int> gamma_model_reconstruction::get_increase_decrease(const gene_family& gf)
{
    clademap<int> size_deltas;
    compute_increase_decrease(_reconstructions[gf.id()].reconstruction, size_deltas);

    return size_deltas;
}

void gamma_model_reconstruction::print_category_likelihoods(std::ostream& ost, const cladevector& order, const std::vector<const gene_family*>& gene_families)
{
    ost << "Family ID\t";
    ostream_iterator<double> lm(ost, "\t");
    copy(_lambda_multipliers.begin(), _lambda_multipliers.end(), lm);
    ost << endl;

    for (auto gf : gene_families)
    {
        ost << gf->id() << '\t';
        auto rc = _reconstructions[gf->id()];
        ostream_iterator<double> ct(ost, "\t");
        copy(rc._category_likelihoods.begin(), rc._category_likelihoods.end(), ct);
        ost << endl;
    }
}

void gamma_model_reconstruction::print_additional_data(const cladevector& order, const std::vector<const gene_family*>& gene_families, std::string output_prefix)
{
    std::ofstream cat_likelihoods(filename("Gamma_category_likelihoods", output_prefix));
    print_category_likelihoods(cat_likelihoods, order, gene_families);

}

int gamma_model_reconstruction::reconstructed_size(const gene_family& family, const clade* clade) const
{
    if (clade->is_leaf())
        return family.get_species_size(clade->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw std::runtime_error("Family " + family.id() + " was not reconstructed");
    auto& c = _reconstructions.at(family.id()).reconstruction;
    if (c.find(clade) == c.end())
        throw std::runtime_error("Clade '" + clade->get_taxon_name() + "' was not reconstructed for family " + family.id());

    return _reconstructions.at(family.id()).reconstruction.at(clade);
}

