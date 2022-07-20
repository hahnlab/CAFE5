#include <cmath>
#include <set>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "doctest.h"
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

struct Report {
    const clade* p_tree = nullptr;
    const lambda* p_lambda = nullptr;
};

std::ostream& operator<<(ostream& ost, const Report& report)
{
    //bool has_pvalues = false;
    //bool has_likelihoods = false;

    ost << "Tree:";
    report.p_tree->write_newick(ost, [](const clade* c)
        {
            ostringstream ost;
            ost << (c->is_leaf() ? c->get_taxon_name() : "") << ":" << c->get_branch_length();
            return ost.str();
        });
    ost << "\n";

    if (report.p_lambda)
    {
        ost << "Lambda:\t";
        auto vals = get_lambda_values(report.p_lambda);
        copy(vals.begin(), vals.end(), ostream_iterator<double>(ost, "\t"));
    }
    ost << "\n";
#if 0
    if (!report.lambda_tree.empty())
    {
        ost << "Lambda tree:\t" << report.lambda_tree << "\n";
    }

    ost << "# IDs of nodes:";

    ost << newick_visualization((pTree)report.aTree) << "\n";

    ost << "# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ";
    for (size_t i = 0; i < report.node_pairs.size(); ++i)
    {
        ost << report.node_pairs[i] << " ";
    }

    ost << "\n";
    ost << "# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (" << report.branch_cutting_output_format[0];

    for (size_t i = 1; i < report.branch_cutting_output_format.size(); i++)
    {
        ost << ", " << report.branch_cutting_output_format[i];
    }
    ost << ")\n";

    write_viterbi(ost, report);

    write_families_header(ost, has_pvalues, has_likelihoods);

    copy(report.family_line_items.begin(), report.family_line_items.end(), ostream_iterator<family_line_item>(ost, "\n"));
#endif
    return ost;
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
                std::ofstream report_file(filename(p_model->name() + "_report", _user_input.output_prefix, "cafe"));
                Report r;
                r.p_tree = data.p_tree;
                r.p_lambda = p_model->get_lambda();
                report_file << r;

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

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("Report writes tree")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Tree:((A:1,B:1):1,(C:1,D:1):1):0");
}

TEST_CASE("Report writes lambdas")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    r.p_lambda = new multiple_lambda(map<std::string, int>(),  { 0.5, 0.3, 0.9 });
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda:\t0.5\t0.3\t0.9");
}
