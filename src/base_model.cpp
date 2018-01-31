#include <cmath>
#include <numeric>
#include <iomanip>

#include "base_model.h"
#include "process.h"
#include "reconstruction_process.h"
#include "matrix_cache.h"

#include "root_equilibrium_distribution.h"

base_model::base_model(lambda* p_lambda, clade *p_tree, vector<gene_family> *p_gene_families,
    int max_family_size, int max_root_family_size, std::map<int, int> * p_rootdist_map, const error_model *p_error_model) :
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_error_model)
{
    if (p_rootdist_map != NULL)
    {
        _rootdist_vec = vectorize_map(p_rootdist_map); // in vector form
        _total_n_families_sim = _rootdist_vec.size();
    }
}

base_model::~base_model()
{
    for (auto proc : processes)
        delete proc;

    for (auto rec_proc : _rec_processes)
        delete rec_proc;
}

simulation_process* base_model::create_simulation_process(int family_number) {
    return new simulation_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec, family_number); // if a single _lambda_multiplier, how do we do it?
}

reconstruction_process* base_model::create_reconstruction_process(int family_number, matrix_cache *p_calc, root_equilibrium_distribution* p_prior) {
    return new reconstruction_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec, 
        &_p_gene_families->at(family_number), p_calc, p_prior);
}

void base_model::start_inference_processes()
{
    processes.clear();
    for (int i = 0; i < _p_gene_families->size(); ++i) {
        inference_process *p_new_process = new inference_process(_ost, _p_lambda, 1.0, _p_tree, _max_family_size, _max_root_family_size, &_p_gene_families->at(i), _rootdist_vec, NULL); // if a single _lambda_multiplier, how do we do it?
        processes.push_back(p_new_process);
    }
}

double base_model::infer_processes(root_equilibrium_distribution *prior) {
#ifdef VERBOSE
    const bool write = true;
#else
    const bool write = false;
#endif
    if (!_p_lambda->is_valid())
    {
        std::cout << "-lnL: " << log(0) << std::endl;
        return -log(0);
    }

    initialize_rootdist_if_necessary();
    prior->initialize(_rootdist_vec);

    results.resize(processes.size());
    std::vector<double> all_families_likelihood(processes.size());

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);

    matrix_cache calc;
    calc.precalculate_matrices(_max_family_size + 1, get_lambda_values(_p_lambda), lengths.result());
    // prune all the families with the same lambda
#pragma omp parallel for

    for (int i = 0; i < processes.size(); ++i) {

        auto partial_likelihood = processes[i]->prune(calc);
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = prior->compute(j);

            double log_full_lk = std::log(partial_likelihood[j]) + std::log(eq_freq);
            full[j] = log_full_lk;

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }

        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)

        results[i] = family_info_stash(i, 0.0, 0.0, 0.0, all_families_likelihood[i], false);
    }

    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    std::cout << "Score (-lnL): " << setw(15) << setprecision(14) << final_likelihood << std::endl;

    return final_likelihood;
}

void base_model::print_results(std::ostream& ost)
{
    ost << "#FamilyID\tLikelihood of Family" << endl;
    for (const auto& r : results)
    {
        ost << r.family_id << "\t" << r.posterior_probability << endl;
    }
}



std::vector<double> base_model::initial_guesses()
{
    branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);
    std::vector<double> result(_p_lambda->count());
    for (auto& i : result)
    {
        i = 1.0 / finder.longest() * unifrnd();
    }
    return result;
}

void base_model::set_current_guesses(double *guesses)
{
    _p_lambda->update(guesses);
    cout << "Lambda: " << *_p_lambda << std::endl;

}

void base_model::reconstruct_ancestral_states(matrix_cache *p_calc, root_equilibrium_distribution* p_prior)
{
    cout << "Starting reconstruction processes for base model" << endl;
    _rec_processes.clear();

    for (int i = 0; i < _p_gene_families->size(); ++i)
    {
        _rec_processes.push_back(create_reconstruction_process(i, p_calc, p_prior));
    }

    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    p_calc->precalculate_matrices(_max_family_size + 1, get_lambda_values(_p_lambda), lengths.result());
    cout << "Base: reconstructing ancestral states - lambda = " << *_p_lambda << endl;
    for (auto p : _rec_processes)
    {
        p->reconstruct();
    }
    cout << "Done!" << endl;
}

void base_model::print_reconstructed_states(std::ostream& ost) {
    auto rec = _rec_processes[0];
    auto order = rec->get_taxa();
    for (auto& it : order) {
        ost << "#" << it->get_taxon_name() << "\n";
    }

    for (int i = 0; i < _rec_processes.size(); ++i) {
        _rec_processes[i]->print_reconstruction(ost, order);
    }
}

