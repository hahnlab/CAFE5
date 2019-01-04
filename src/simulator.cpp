#include "simulator.h"
#include "user_data.h"
#include "core.h"
#include "matrix_cache.h"

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

trial* simulator::create_trial(model *p_model, const root_distribution& rd, int family_number, const user_data& data, const matrix_cache& cache) {

    if (data.p_lambda == NULL)
        throw std::runtime_error("No lambda specified for simulation");

    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    int max_family_size_sim;

    auto *result = new trial();

    if (data.rootdist.empty()) {
        max_family_size_sim = 100;
    }
    else {
        max_family_size_sim = 2 * rd.max();
    }

    (*result)[data.p_tree] = select_root_size(data, rd, family_number);

    unique_ptr<lambda> sim_lambda(p_model->get_simulation_lambda(data));
    random_familysize_setter rfs(result, max_family_size_sim, sim_lambda.get(), data.p_error_model, cache);

    data.p_tree->apply_prefix_order(rfs); // this is where the () overload of random_familysize_setter is used

    return result;
}


void simulator::simulate_processes(model *p_model, std::vector<trial *>& results) {

    int max_size = p_model->get_max_simulation_size();
    size_t rootdist_sz = p_model->get_rootdist_size();
    if (rootdist_sz > 0)
    {
        results.resize(rootdist_sz);
    }
    else
    {
        results.resize(_user_input.nsims);
    }

#ifndef SILENT
    cout << "Simulating " << results.size() << " families for model " << p_model->name() << endl;
#endif

    root_distribution rd;

    if (data.rootdist.empty()) {
        rd.vectorize_increasing(100);
    }
    else {
        rd.vectorize(data.rootdist);
    }

    for (size_t i = 0; i < results.size(); i+=50)
    {
        p_model->perturb_lambda();

        matrix_cache cache(max_size);
        p_model->prepare_matrices_for_simulation(cache);

#ifndef SILENT
        cout << "Matrices complete\n";
        cache.warn_on_saturation(cerr);
#endif
        int n = 0;

        auto end_it = i + 50 > results.size() ? results.end() : results.begin() + i + 50;
        generate(results.begin()+i, end_it, [this, p_model, i, &rd, &cache, &n]() mutable {
            return create_trial(p_model, rd, i+n++, data, cache);
        });
    }
}


/// Simulate
/// \callgraph
void simulator::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    cout << "Simulating with " << models.size() << " model(s)" << endl;

    std::vector<const clade *> order;
    auto fn = [&order](const clade *c) { order.push_back(c); };
    data.p_tree->apply_reverse_level_order(fn);

    for (auto p_model : models) {

        std::vector<trial *> results;

        simulate_processes(p_model, results);

        string truth_fname = filename("simulation_truth", my_input_parameters.output_prefix);
        std::ofstream ofst(truth_fname);
        cout << "Writing to " << truth_fname << endl;
        print_simulations(ofst, true, results);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        cout << "Writing to " << fname << endl;
        print_simulations(ofst2, false, results);
    }
}


void simulator::print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<trial *>& results) {

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
