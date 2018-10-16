#include <cmath>
#include <set>

#include "utils.h"
#include "io.h"
#include "execute.h"
#include "probability.h"
#include "lambda.h"
#include "core.h"
#include "gamma_core.h"
#include "root_equilibrium_distribution.h"
#include "chisquare.h"
#include "matrix_cache.h"
#include "base_model.h"
#include "user_data.h"

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


std::string filename(std::string base, std::string suffix)
{
    return base + (suffix.empty() ? "" : "_") + suffix + ".txt";
}

void estimator::compute(std::vector<model *>& models, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size)
{
    std::ofstream results_file(filename("results", my_input_parameters.output_prefix));

    std::ofstream likelihoods_file(filename("family_lks", my_input_parameters.output_prefix));

    std::vector<double> model_likelihoods(models.size());
    for (int i = 0; i < models.size(); ++i) {
        cout << endl << "Starting inference processes for " << models[i]->name() << " model" << endl;
        models[i]->start_inference_processes();

        cout << endl << "Inferring processes for " << models[i]->name() << " model" << endl;
        double result = models[i]->infer_processes(data.p_prior.get());
        models[i]->write_vital_statistics(results_file, result);

        models[i]->write_family_likelihoods(likelihoods_file);

        model_likelihoods[i] = result;
    }

    if (model_likelihoods.size() == 2)
    {
        cout << "PValue = " << (1.0 - chi2cdf(2 * (model_likelihoods[1] - model_likelihoods[0]), 1.0));
    }
}

bool compare_result(const optimizer::result& a, const optimizer::result& b)
{
    return a.score < b.score;
}

void estimator::estimate_missing_variables(std::vector<model *>& models, user_data& data)
{
    if (data.p_tree == NULL)
    {
        throw runtime_error("No tree specified for lambda estimation");
    }
    for (model* p_model : models) {

        unique_ptr<optimizer_scorer> scorer(p_model->get_lambda_optimizer(data));
        if (scorer.get() == nullptr)
            continue;   // nothing to be optimized

        optimizer opt(scorer.get());

        if (data.p_error_model)
        {
            // try several different initialization points and optimize them
            vector<optimizer::result> results;
            for (double epsilon = 0.05; epsilon < .5; epsilon += .1)
            {
                data.p_error_model->update_single_epsilon(epsilon);
                results.push_back(opt.optimize());
            }
            auto best = min_element(results.begin(), results.end(), compare_result);
            scorer->finalize(&best->values[0]);
            cout << "Final score: " << best->score << ", Lambda: " << best->values[0] << ", Epsilon: " << best->values[1] * 2 << endl;
        }
        else
        {
            auto result = opt.optimize();
            scorer->finalize(&result.values[0]);
        }
    }
    if (data.p_lambda == nullptr)
        data.p_lambda = models[0]->get_lambda();

}

void estimator::execute(std::vector<model *>& models)
{
    estimate_missing_variables(models, data);

    compute(models, _user_input, data.max_family_size, data.max_root_family_size);

    auto pvalues = compute_pvalues(data, 1000);

    matrix_cache cache(data.max_family_size + 1);
    for (model* p_model : models) {
        std::unique_ptr<reconstruction> rec(p_model->reconstruct_ancestral_states(&cache, data.p_prior.get()));
        rec->write_results(p_model, _user_input.output_prefix, pvalues);
    }
}

double pvalue(double v, const vector<double>& conddist)
{
    int idx = conddist.size() - 1;

    auto bound = std::upper_bound(conddist.begin(), conddist.end(), v);
    if (bound != conddist.end())
    {
        idx = bound - conddist.begin();
    }
    return  idx / (double)conddist.size();
}

class pvalue_calculator
{
    int _max_family_size;
    int _max_root_family_size;
    const lambda *_p_lambda; 
    const clade* _p_tree;
    matrix_cache *_p_matrix_cache;
    const std::vector<std::vector<double> >* _p_conditional_distribution;
public:
    pvalue_calculator(int max_family_size, int max_root_family_size, const lambda *p_lambda, const clade* p_tree, matrix_cache *p_matrix_cache, const std::vector<std::vector<double> >* p_conditional_distribution) :
        _max_family_size(max_family_size), 
        _max_root_family_size(max_root_family_size),
        _p_lambda(p_lambda),
        _p_tree(p_tree),
        _p_matrix_cache(p_matrix_cache),
        _p_conditional_distribution(p_conditional_distribution)
    {

    }

    double operator()(const gene_family& gf);
};

double pvalue_calculator::operator()(const gene_family& gf)
{
    int max = gf.get_max_size();

    likelihood_computer pruner(_max_root_family_size, _max_family_size, _p_lambda, gf, *_p_matrix_cache, NULL);
    _p_tree->apply_reverse_level_order(pruner);
    auto lh = pruner.get_likelihoods(_p_tree);

    double observed_max_likelihood = pruner.max_likelihood(_p_tree);	// max value but do we need a posteriori value instead?

    vector<double> pvalues(_max_root_family_size);
    for (int s = 0; s < _max_root_family_size; s++)
    {
        pvalues[s] = pvalue(observed_max_likelihood, _p_conditional_distribution->at(s));
    }
    return *max_element(pvalues.begin(), pvalues.end());
}

// find_fast_families under base model through simulations (if we reject gamma)
vector<double> estimator::compute_pvalues(const user_data& data, int number_of_simulations)
{
    cout << "Computing pvalues..." << flush;

    matrix_cache cache(max(data.max_family_size, data.max_root_family_size) + 1);
    branch_length_finder lengths;
    data.p_tree->apply_prefix_order(lengths);
    cache.precalculate_matrices(get_lambda_values(data.p_lambda), lengths.result());

    auto cd = get_conditional_distribution_matrix(data.p_tree, data.max_root_family_size, data.max_family_size, number_of_simulations, data.p_lambda, cache);

    vector<double> result(data.gene_families.size());

    pvalue_calculator calculator(data.max_family_size, data.max_root_family_size, data.p_lambda, data.p_tree, &cache, &cd);
    transform(data.gene_families.begin(), data.gene_families.end(), result.begin(), calculator);

    cout << "done!\n";

    return result;
}

void chisquare_compare::execute(std::vector<model *>&)
{
    vector<string> chistrings = tokenize_str(_values, ',');
    vector<double> chis(chistrings.size());

    // transform is like R's apply (vector lambdas takes the outputs, here we are making doubles from strings
    transform(chistrings.begin(), chistrings.end(), chis.begin(),
        [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
    );

    double degrees_of_freedom = chis[2];
    cout << "PValue = " << 1.0 - chi2cdf(2 * (chis[1] - chis[0]), degrees_of_freedom) << std::endl;
}

void simulator::execute(std::vector<model *>& models)
{
    if (models.empty())
        throw std::runtime_error("Not enough information to specify a model");

    if (_user_input.rootdist.empty())
    {
        throw std::runtime_error("No root distribution specified"); // cannot simulate without a rootdist (TODO)
    }

    // -s is provided an argument (-f is not), using -i to obtain root eq freq distr'n
    if (_user_input.nsims != 0) {
        throw std::runtime_error("A specified number of simulations with a root distribution is not supported");
        // place holder for estimating poisson lambda if -p, or using uniform as root eq freq distr'n
    }

    // -f is provided (-s does not have an argument), not using -i
    else if (!_user_input.rootdist.empty()) {
        simulate(models, _user_input);
    }
}

/// Simulate
/// \callgraph
void simulator::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    cout << "Simulating with " << models.size() << " model(s)" << endl;

    for (int i = 0; i < models.size(); ++i) {

        cout << "Simulating " << models[i]->get_total_n_families_sim() << " families for model " << models[i]->name() << endl;

        // lambda_multipliers and lambda_bins will not be harcoded in the future
        if (my_input_parameters.rootdist.empty())
            models[i]->set_total_n_families_sim(my_input_parameters.nsims);

        models[i]->start_sim_processes();

        models[i]->simulate_processes();

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst(fname);
        cout << "Writing to " << fname << endl;
        models[i]->print_processes(ofst);
    }
}

