#include <cmath>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "base_model.h"
#include "gene_family_reconstructor.h"
#include "matrix_cache.h"
#include "gene_family.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "optimizer_scorer.h"
#include "simulator.h"

#if defined __INTEL_COMPILER
#include <pstl/execution>
#include <pstl/algorithm> 
#elif defined __PGI
#include <pstl/execution>
#include <pstl/algorithm> 
#elif defined __llvm__
#include <pstl/execution>
#include <pstl/algorithm> 
#elif defined _CRAYC
#include <pstl/execution>
#include <pstl/algorithm> 
#elif defined __GNUC__
#include <execution>
#endif

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

double base_model::infer_family_likelihoods(const root_equilibrium_distribution &prior, const lambda *p_lambda) {
    _monitor.Event_InferenceAttempt_Started();

    if (!_p_lambda->is_valid())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    results.resize(_p_gene_families->size());
    std::vector<double> all_families_likelihood(_p_gene_families->size());

    matrix_cache calc(max(_max_root_family_size, _max_family_size) + 1);
    calc.precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());

    vector<vector<double>> partial_likelihoods(_p_gene_families->size());
#ifdef USE_STDLIB_PARALLEL
    par_timer.start("stdlib: inference prune (Base)");
    transform(std::execution::par, _p_gene_families->begin(), _p_gene_families->end(), partial_likelihoods.begin(), 
		[&](const gene_family &f) {
				return inference_prune(f, calc, _p_lambda, _p_error_model, _p_tree, 1.0, _max_root_family_size, _max_family_size);
			});
    par_timer.stop("stdlib: inference prune (Base)");
#else
    par_timer.start("OMP: inference prune (Base)");
#pragma omp parallel for
        for (size_t i = 0; i < _p_gene_families->size(); ++i) {
            if (references[i] == i)
                partial_likelihoods[i] = inference_prune(_p_gene_families->at(i), calc, _p_lambda, _p_error_model, _p_tree, 1.0, _max_root_family_size, _max_family_size);
                // probabilities of various family sizes
    }
    par_timer.stop("OMP: inference prune (Base)");
#endif

    // prune all the families with the same lambda
#ifdef USE_STDLIB_PARALLEL
    par_timer.start("stdlib: compute likelihood (Base)");
    std::map<string, int> refmap;
    for (int i = 0; i < _p_gene_families->size(); ++i)
    {
        refmap[_p_gene_families->at(i).id()] = references[i];
    }
    transform(std::execution::par, _p_gene_families->begin(), _p_gene_families->end(), results.begin(), [&](const gene_family& gf) {

        auto& partial_likelihood = partial_likelihoods[refmap[gf.id()]];
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = prior.compute(j);

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }

        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        double mx = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
                                                                             // cout << i << " contribution " << scientific << all_families_likelihood[i] << endl;

        return family_info_stash(gf.id(), 0.0, 0.0, 0.0, mx, false);
        });
    par_timer.stop("stdlib: compute likelihood (Base)");
#else
    par_timer.start("OMP: compute likelihood (Base)");
#pragma omp parallel for
    for (size_t i = 0; i < _p_gene_families->size(); ++i) {

        auto& partial_likelihood = partial_likelihoods[references[i]];
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = prior.compute(j);

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }

        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
                                                                             // cout << i << " contribution " << scientific << all_families_likelihood[i] << endl;

        results[i] = family_info_stash(_p_gene_families->at(i).id(), 0.0, 0.0, 0.0, all_families_likelihood[i], false);
    }
    par_timer.stop("OMP: compute likelihood (Base)");
#endif

    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

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

inference_optimizer_scorer *base_model::get_lambda_optimizer(const user_data& data)
{
    if (data.p_lambda != NULL)  // already have a lambda, nothing we want to optimize
        return nullptr;

    initialize_lambda(data.p_lambda_tree);

    auto lengths = _p_tree->get_branch_lengths();
    auto longest_branch = *max_element(lengths.begin(), lengths.end());

    if (_p_error_model && !data.p_error_model)
    {
        return new lambda_epsilon_optimizer(this, _p_error_model, &data.prior, data.rootdist, _p_lambda, longest_branch);
    }
    else
    {
        return new lambda_optimizer(_p_lambda, this, &data.prior, longest_branch);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(const vector<gene_family>& families, matrix_cache *p_calc, root_equilibrium_distribution* p_prior)
{
    LOG(INFO) << "Starting reconstruction processes for Base model";

    auto result = new base_model_reconstruction();

    p_calc->precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());

    pupko_reconstructor::pupko_data data(families.size(), _p_tree, _max_family_size, _max_root_family_size);

    for (size_t i = 0; i < families.size(); ++i)
    {
        clademap<int> &rc = result->_reconstructions[families[i].id()];
        _p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c] = 0;
            });
    }

#pragma omp parallel for
    for (size_t i = 0; i< families.size(); ++i)
    {
        pupko_reconstructor::reconstruct_gene_family(_p_lambda, _p_tree, &families[i], p_calc, p_prior, result->_reconstructions[families[i].id()], data.C(i), data.L(i));
    }

    size_t success = count_if(data.v_all_node_Ls.begin(), data.v_all_node_Ls.end(), [this](const clademap<std::vector<double>>& L) 
        { 
            return *max_element( L.at(_p_tree).begin(), L.at(_p_tree).end()) > 0; 
        });

    if (success != families.size())
    {
        LOG(WARNING) << "Failed to reconstruct " << families.size() - success << " families" << endl;
    }

    LOG(INFO) << "Done!\n";

    return result;
}

void base_model::prepare_matrices_for_simulation(matrix_cache& cache)
{
    unique_ptr<lambda> perturbed_lambda(get_simulation_lambda());
    cache.precalculate_matrices(get_lambda_values(_p_lambda), _p_tree->get_branch_lengths());
}

lambda* base_model::get_simulation_lambda()
{
    return _p_lambda->multiply(simulation_lambda_multiplier);
}

int base_model_reconstruction::get_node_count(const gene_family& family, const clade *c) const
{
    if (c->is_leaf())
        return family.get_species_size(c->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw runtime_error("Family " + family.id() + " not found in reconstruction");

    return _reconstructions.at(family.id()).at(c);
}
