#include <cmath>
#include <set>

#include "execute.h"
#include "core.h"
#include "user_data.h"
#include "chisquare.h"
#include "optimizer_scorer.h"
#include "matrix_cache.h"
#include "optimizer.h"

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


void estimator::compute(std::vector<model *>& models, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size)
{
    std::ofstream results_file(filename("results", my_input_parameters.output_prefix));

    std::ofstream likelihoods_file(filename("family_lks", my_input_parameters.output_prefix));

    std::vector<double> model_likelihoods(models.size());
    for (size_t i = 0; i < models.size(); ++i) {
        cout << endl << "Starting inference processes for " << models[i]->name() << " model" << endl;
        models[i]->start_inference_processes(models[i]->get_lambda());

        cout << endl << "Inferring processes for " << models[i]->name() << " model" << endl;
        double result = models[i]->infer_processes(data.p_prior.get(), data.rootdist);
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

        unique_ptr<inference_optimizer_scorer> scorer(p_model->get_lambda_optimizer(data));
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

void estimator::estimate_lambda_per_family(model *p_model, ostream& ost)
{
    vector<lambda*> result(data.gene_families.size());
    std::transform(data.gene_families.begin(), data.gene_families.end(), result.begin(),
        [this, p_model](gene_family& fam)
    {
#ifndef SILENT
        cout << "Estimating for " << fam.id() << endl;
#endif
        vector<gene_family> v({ fam });
        vector<model *> models{ p_model };
        p_model->set_families(&v);
        data.p_lambda = nullptr;
        estimate_missing_variables(models, data);
        return p_model->get_lambda();
    });
    std::transform(data.gene_families.begin(), data.gene_families.end(), result.begin(),
        ostream_iterator<string>(ost, "\n"),
        [](const gene_family& fam, lambda* lambda)
    {
        return fam.id() + '\t' + lambda->to_string();
    });
    for (auto r : result) delete r; // TODO: use unique_ptrs in result for exception safety

}

void estimator::execute(std::vector<model *>& models)
{
    if (_user_input.lambda_per_family)
    {
        auto p_model = models[0];   // no support for multiple models
        std::ofstream results_file(filename(p_model->name() + "_lambda_per_family", _user_input.output_prefix));
        estimate_lambda_per_family(p_model, results_file);
    }
    else
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

// find_fast_families under base model through simulations (if we reject gamma)
vector<double> estimator::compute_pvalues(const user_data& data, int number_of_simulations) const
{
#ifndef SILENT
    cout << "Computing pvalues..." << flush;
#endif

    matrix_cache cache(max(data.max_family_size, data.max_root_family_size) + 1);
    branch_length_finder lengths;
    data.p_tree->apply_prefix_order(lengths);
    cache.precalculate_matrices(get_lambda_values(data.p_lambda), lengths.result());

    std::vector<std::vector<double> > cd(data.max_root_family_size);
    for (int i = 0; i < data.max_root_family_size; ++i)
    {
        cd[i] = get_random_probabilities(data.p_tree, number_of_simulations, i, data.max_family_size, data.max_root_family_size, data.p_lambda, cache, NULL);
    }

    vector<double> result(data.gene_families.size());

    transform(data.gene_families.begin(), data.gene_families.end(), result.begin(), [&data, &cache, &cd](const gene_family& gf)
    {
        likelihood_computer pruner(data.max_root_family_size, data.max_family_size, data.p_lambda, gf, cache, NULL);
        pruner.initialize_memory(data.p_tree);
        data.p_tree->apply_reverse_level_order(pruner);
        auto lh = pruner.get_likelihoods(data.p_tree);

        double observed_max_likelihood = pruner.max_likelihood(data.p_tree);	// max value but do we need a posteriori value instead?

        vector<double> pvalues(data.max_root_family_size);
        for (int s = 0; s < data.max_root_family_size; s++)
        {
            pvalues[s] = pvalue(observed_max_likelihood, cd[s]);
        }
        return *max_element(pvalues.begin(), pvalues.end());
    });

#ifndef SILENT
    cout << "done!\n";
#endif
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

