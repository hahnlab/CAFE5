#include <numeric>
#include <algorithm>
#include <fstream>

#include "simulator.h"
#include "user_data.h"
#include "core.h"
#include "matrix_cache.h"

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

int select_root_size(const user_data& data, const root_distribution& rd, int family_number)
{
    if (data.rootdist.empty()) {
        return rd.select_randomly(); // getting a random root size from the provided (core's) root distribution
    }
    else {
        return rd.at(family_number);
    }
}

clademap<int>* simulator::create_trial(const lambda *p_lambda, const root_distribution& rd, int family_number, const matrix_cache& cache) {

    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    int max_family_size_sim;

    auto *result = new clademap<int>();

    if (data.rootdist.empty()) {
        max_family_size_sim = 100;
    }
    else {
        max_family_size_sim = 2 * rd.max();
    }

    (*result)[data.p_tree] = select_root_size(data, rd, family_number);


    auto fn = [&](const clade *c)
    {
        set_weighted_random_family_size(c, result, p_lambda, data.p_error_model, max_family_size_sim, cache);
    };

    data.p_tree->apply_prefix_order(fn);

    return result;
}


void simulator::simulate_processes(model *p_model, std::vector<clademap<int> *>& results) {

    root_distribution rd;
    int max_size;
    if (data.rootdist.empty())
    {
        results.resize(_user_input.nsims);
        max_size = 100;
        rd.vectorize_increasing(max_size);
    }
    else
    {
        rd.vectorize(data.rootdist);
        if (_user_input.nsims > 0)
        {
            rd.pare(_user_input.nsims);
        }
        results.resize(rd.size());
        max_size = 2 * rd.max();
    }

    if (!quiet)
        cout << endl << "Simulating " << results.size() << " families for model " << p_model->name() << endl << endl;

    for (size_t i = 0; i < results.size(); i+= LAMBDA_PERTURBATION_STEP_SIZE)
    {
        p_model->perturb_lambda();
        unique_ptr<lambda> sim_lambda(p_model->get_simulation_lambda());
        
        matrix_cache cache(max_size);
        //cache.precalculate_matrices(get_lambda_values(sim_lambda.get()), this->data.p_tree->get_branch_lengths());
        p_model->prepare_matrices_for_simulation(cache);

        if (!quiet)
            cache.warn_on_saturation(cerr);

        int n = 0;

        auto end_it = i + LAMBDA_PERTURBATION_STEP_SIZE > results.size() ? results.end() : results.begin() + i + LAMBDA_PERTURBATION_STEP_SIZE;
        generate(results.begin()+i, end_it, [this, &sim_lambda, i, &rd, &cache, &n]() mutable {
            return create_trial(sim_lambda.get(), rd, i+n++, cache);
        });
    }
}

extern void write_average_multiplier(std::ostream& ost);

/// Simulate
/// \callgraph
void simulator::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    cout << endl << "Simulating with " << models.size() << " model(s)" << endl; 

	if (data.p_tree == nullptr)
		throw std::runtime_error("No tree specified for simulations");

    std::vector<const clade *> order;
    data.p_tree->apply_reverse_level_order([&order](const clade* c) { order.push_back(c); });

    string dir = my_input_parameters.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);

    for (auto p_model : models) {

        std::vector<clademap<int> *> results;

        simulate_processes(p_model, results);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        print_simulations(ofst2, false, results);
        if (!quiet)
            cout << "Simulated values written to " << fname << endl;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
        print_simulations(ofst, true, results);
        if (!quiet)
            cout << "Simulated values (including internal nodes) written to " << truth_fname << endl;

        if (my_input_parameters.fixed_lambda > 0)
        {
            write_average_multiplier(cout);
        }

    }
}


void simulator::print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<clademap<int> *>& results) {

    std::vector<const clade *> order;
    auto fn = [&order](const clade *c) { order.push_back(c); };
    data.p_tree->apply_reverse_level_order(fn);

    if (results.empty())
    {
        cerr << "No simulations created" << endl;
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
        auto& fam = *results[j];
        // Printing gene counts
        ost << "NULL\tsimfam" << j;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (order[i]->is_leaf() || include_internal_nodes)
            {
                ost << '\t';
                ost << fam[order[i]];
            }
        }
        ost << endl;
    }
}
