#include <numeric>
#include <algorithm>
#include <fstream>
#include <random>

#include "easylogging++.h"

#include "simulator.h"
#include "user_data.h"
#include "core.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"

extern std::mt19937 randomizer_engine; // seeding random number engine

simulator::simulator(user_data& d, const input_parameters& ui) : action(d, ui)
{
#ifdef SILENT
    quiet = true;
#endif
}

void simulator::execute(std::vector<model *>& models)
{
    simulate(models, _user_input);
}

simulated_family simulator::create_trial(const lambda *p_lambda, int family_number, const matrix_cache& cache) {

    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    simulated_family result;
    result.lambda = get_lambda_values(p_lambda)[0];

    int i = 0;
    for (i = 0; i<50; ++i)
    {
        result.values[data.p_tree] = data.prior.select_root_size(family_number);

        data.p_tree->apply_prefix_order([&](const clade* c)
            {
                set_weighted_random_family_size(c, &result.values, p_lambda, data.p_error_model, data.max_family_size, cache);
            });

        gene_family gf;
        gf.init_from_clademap(result.values);
        if (gf.exists_at_root(data.p_tree))
            break;
    }
    if (i >= 50)
    {
        LOG(WARNING) << "Failed to create a family that would exist at the root\n";
    }

    return result;
}

void simulator::simulate_processes(model *p_model, std::vector<simulated_family>& results) {

    if (_user_input.nsims > 0)
    {
        results.resize(_user_input.nsims);
    }
    else
    {
        results.resize(accumulate(data.rootdist.begin(), data.rootdist.end(), 0,
            [](int acc, std::pair<int, int> p) { return (acc + p.second); }));
    }

    LOG(INFO) << "Simulating " << results.size() << " families for model " << p_model->name();

    for (size_t i = 0; i < results.size(); i+= LAMBDA_PERTURBATION_STEP_SIZE)
    {
        p_model->perturb_lambda();
        unique_ptr<lambda> sim_lambda(p_model->get_simulation_lambda());
        
        matrix_cache cache(data.max_root_family_size+1);
        //cache.precalculate_matrices(get_lambda_values(sim_lambda.get()), this->data.p_tree->get_branch_lengths());
        p_model->prepare_matrices_for_simulation(cache);

        if (!quiet)
            cache.warn_on_saturation(cerr);

        int n = 0;

        auto end_it = i + LAMBDA_PERTURBATION_STEP_SIZE > results.size() ? results.end() : results.begin() + i + LAMBDA_PERTURBATION_STEP_SIZE;
        generate(results.begin()+i, end_it, [this, &sim_lambda, i, &cache, &n]() mutable {
            return create_trial(sim_lambda.get(), i+n++, cache);
        });
    }
}

extern void write_average_multiplier(std::ostream& ost);

/// Simulate
/// \callgraph
void simulator::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    LOG(INFO) << "Simulating with " << models.size() << " model(s)";

	if (data.p_tree == nullptr)
		throw std::runtime_error("No tree specified for simulations");

    std::vector<const clade *> order;
    data.p_tree->apply_reverse_level_order([&order](const clade* c) { order.push_back(c); });

    string dir = my_input_parameters.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);

    for (auto p_model : models) {

        std::vector<simulated_family> results;

        simulate_processes(p_model, results);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        print_simulations(ofst2, false, results);
        LOG(INFO) << "Simulated values written to " << fname << endl;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
        print_simulations(ofst, true, results);
        LOG(INFO) << "Simulated values (including internal nodes) written to " << truth_fname << endl;

        if (my_input_parameters.fixed_lambda > 0)
        {
            write_average_multiplier(cout);
        }

    }
}


void simulator::print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<simulated_family>& results) {

    std::vector<const clade *> order;
    auto fn = [&order](const clade *c) { order.push_back(c); };
    data.p_tree->apply_reverse_level_order(fn);

    if (results.empty())
    {
        LOG(ERROR) << "No simulations created" << endl;
        return;
    }
    ost << "DESC\tFID";
    for (size_t i = 0; i < order.size(); ++i)
    {
        if (order[i]->is_leaf())
            ost << '\t' << order[i]->get_taxon_name();
        else if (include_internal_nodes)
            ost << '\t' << i;

    }
    ost << endl;

    for (size_t j = 0; j < results.size(); ++j) {
        auto& fam = results[j];
        // Printing gene counts
        ost << "L" << fam.lambda << "\tsimfam" << j;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (order[i]->is_leaf() || include_internal_nodes)
            {
                ost << '\t';
                ost << fam.values.at(order[i]);
            }
        }
        ost << endl;
    }
}
