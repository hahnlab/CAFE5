#include <cmath>
#include <numeric>
#include <limits>
#include <random>

#include "base_model.h"
#include "process.h"
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

base_model::~base_model()
{
    for (auto proc : processes)
        delete proc;

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

void base_model::start_inference_processes(lambda *p_lambda)
{
    for (auto proc : processes)
        delete proc;
    processes.clear();

    processes.resize(_p_gene_families->size());
    for (size_t i = 0; i <_p_gene_families->size(); ++i) {
        if (references[i] == i)
            processes[i] = new inference_process(_ost, p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _p_error_model); // if a single _lambda_multiplier, how do we do it?
    }
}

double base_model::infer_processes(root_equilibrium_distribution *prior, const std::map<int, int>& root_distribution_map) {
    if (!_p_lambda->is_valid())
    {
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

    results.resize(processes.size());
    std::vector<double> all_families_likelihood(processes.size());

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);

    matrix_cache calc(max(_max_root_family_size, _max_family_size) + 1);
    calc.precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());

    vector<vector<double>> partial_likelihoods(processes.size());
#pragma omp parallel for
    for (size_t i = 0; i < processes.size(); ++i) {
        if (processes[i])
            partial_likelihoods[i] = processes[i]->prune(calc);    // probabilities of various family sizes
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
#ifndef SILENT
    cout << "Starting reconstruction processes for base model" << endl;
#endif

    base_model_reconstruction *result = new base_model_reconstruction();

    result->_rec_processes.resize(_p_gene_families->size());
    transform(_p_gene_families->begin(), _p_gene_families->end(), result->_rec_processes.begin(), [this, p_calc, p_prior](const gene_family &family)
    {
        return new gene_family_reconstructor(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size,
            &family, p_calc, p_prior);
    });


    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    p_calc->precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());

#ifndef SILENT
    cout << "Base: reconstructing ancestral states - lambda = " << *_p_lambda << endl;
#endif

    for (auto p : result->_rec_processes)
    {
        p->reconstruct();
    }
#ifndef SILENT
    cout << "Done!" << endl;
#endif

    return result;
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

base_model_reconstruction::~base_model_reconstruction()
{
    for (auto rec_proc : _rec_processes)
        delete rec_proc;
}

void base_model_reconstruction::print_reconstructed_states(std::ostream& ost) {
    if (_rec_processes.empty())
        return;

    auto rec = _rec_processes[0];
    auto order = rec->get_taxa();

    ost << "#NEXUS\nBEGIN TREES;\n";
    for (auto item : _rec_processes)
    {
        item->print_reconstruction(ost, order);
    }
    ost << "END;\n";
}

void base_model_reconstruction::print_increases_decreases_by_family(std::ostream& ost, const std::vector<double>& pvalues) {
    ::print_increases_decreases_by_family(ost, _rec_processes, pvalues);
}

void base_model_reconstruction::print_increases_decreases_by_clade(std::ostream& ost) {
    ::print_increases_decreases_by_clade(ost, _rec_processes);
}

