#include <cmath>
#include <set>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "easylogging++.h"

#include "execute.h"
#include "core.h"
#include "user_data.h"
#include "chisquare.h"
#include "optimizer_scorer.h"
#include "matrix_cache.h"
#include "optimizer.h"
#include "gene_family_reconstructor.h"
#include "error_model.h"
#include "likelihood_ratio.h"
#include "io.h"
#include "root_equilibrium_distribution.h"

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


void estimator::write_error_model_if_specified(const input_parameters& my_input_parameters, const model * p_model)
{
    if (my_input_parameters.use_error_model)
    {
        ofstream errmodel(filename(p_model->name() + "_error_model", _user_input.output_prefix));
        if (data.p_error_model)
        {
            /// user specified an error model, write that out to the results directory
            write_error_model_file(errmodel, *data.p_error_model);
        }
        else
        {
            /// user did not specify an error model, write out the estimated one or a default
            p_model->write_error_model(errmodel);
        }
    }
}

void estimator::compute(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    std::vector<double> model_likelihoods(models.size());
    for (size_t i = 0; i < models.size(); ++i) {
        LOG(INFO) << "Inferring processes for " << models[i]->name() << " model";

        double result = models[i]->infer_family_likelihoods(data.prior, models[i]->get_lambda());
        std::ofstream results_file(filename(models[i]->name() + "_results", my_input_parameters.output_prefix));
        models[i]->write_vital_statistics(results_file, result);

        std::ofstream likelihoods_file(filename(models[i]->name() + "_family_likelihoods", my_input_parameters.output_prefix));
        models[i]->write_family_likelihoods(likelihoods_file);

        write_error_model_if_specified(my_input_parameters, models[i]);

        model_likelihoods[i] = result;
    }

    auto lengths = data.p_tree->get_branch_lengths();
    auto longest_branch = *max_element(lengths.begin(), lengths.end());
    auto max_lambda = 1 / longest_branch;

    LOG(INFO) << "Maximum possible lambda for this topology: " << max_lambda;

    if (model_likelihoods.size() == 2)
    {
        LOG(INFO) << "PValue = " << (1.0 - chi2cdf(2 * (model_likelihoods[1] - model_likelihoods[0]), 1.0));
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

        auto result = opt.optimize(_user_input.optimizer_params);
        scorer->finalize(&result.values[0]);

        LOG(INFO) << p_model->get_monitor();
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
            LOG(INFO) << "Estimating for " << fam.id() << endl;
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

/*! Calls estimate_lambda_per_family if the user has set that parameter, otherwise
    calls \ref estimate_missing_variables; \ref compute; \ref compute_pvalues, and \ref model::reconstruct_ancestral_states */
void estimator::execute(std::vector<model *>& models)
{
    string dir = _user_input.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);
    
    if (_user_input.lambda_per_family)
    {
        auto p_model = models[0];   // no support for multiple models
        std::ofstream results_file(filename(p_model->name() + "_lambda_per_family", _user_input.output_prefix));
        estimate_lambda_per_family(p_model, results_file);
    }
    else
    {
        try
        {
            estimate_missing_variables(models, data);

            compute(models, _user_input);

            matrix_cache cache(data.max_family_size + 1);
            for (model* p_model : models) {

                /// For Gamma models, we tried using the most rapidly changing lambda multiplier here, but that
                /// caused issues in the pvalue calculation. It should be best to use the original lambda
                /// instead
                matrix_cache cache(max(data.max_family_size, data.max_root_family_size) + 100);
                cache.precalculate_matrices(get_lambda_values(p_model->get_lambda()), data.p_tree->get_branch_lengths());

                pvalue_parameters p = { data.p_tree, p_model->get_lambda(), data.max_family_size, data.max_root_family_size, cache };
                auto pvalues = compute_pvalues(p, data.gene_families, 1000 );

                std::unique_ptr<reconstruction> rec(p_model->reconstruct_ancestral_states(data.gene_families, &cache, &data.prior));

                branch_probabilities probs;

                for (size_t i = 0; i<data.gene_families.size(); ++i)
                {
                    if (pvalues[i] < _user_input.pvalue)
                    {
                        for_each(data.p_tree->reverse_level_begin(), data.p_tree->reverse_level_end(), [&](const clade* c) {
                            probs.set(data.gene_families[i], c, compute_viterbi_sum(c, data.gene_families[i], rec.get(), data.max_family_size, cache, p_model->get_lambda()));
                            });
                    }
                }

#ifdef RUN_LHRTEST
                LikelihoodRatioTest::lhr_for_diff_lambdas(data, p_model);
#endif
                rec->write_results(p_model->name(), _user_input.output_prefix, data.p_tree, data.gene_families, pvalues, _user_input.pvalue, probs);
            }
        }
        catch (const OptimizerInitializationFailure& e )
        {
            initialization_failure_advice(cerr, data.gene_families);
            throw;
        }
    }
}

//Calculate the difference between the Max and Min count for each family, report the 20 families with the largest difference.
void initialization_failure_advice(std::ostream& ost, const std::vector<gene_family>& families)
{
    std::vector<std::pair<std::string, int>> m;
    transform(families.begin(), families.end(), std::inserter(m, m.end()),
        [](const gene_family& gf) { return std::make_pair(gf.id(), gf.species_size_differential()); });
    auto compare = [](const std::pair<string, int>& a, const std::pair<string, int>& b) { return a.second > b.second; };
    sort(m.begin(), m.end(), compare);
    if (m.size() > 20)
        m.resize(20);

    ost << "\nFamilies with largest size differentials:\n";
    for (auto& t : m)
        ost << t.first << ": " << t.second << "\n";
    ost << "\nYou may want to try removing the top few families with the largest difference\nbetween the max and min counts and then re-run the analysis.\n\n";
}
