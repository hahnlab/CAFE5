#include "core.h"
#include "matrix_cache.h"
#include "user_data.h"
#include "chisquare.h"
#include "optimizer_scorer.h"

#include <memory>
#include <cmath>
#include <algorithm>

namespace LikelihoodRatioTest
{
    using namespace std;

    clade * update_branchlength(const clade * p_tree, double bl_augment, int t)
    {
        return new clade(*p_tree, nullptr, [&](const clade& c) {
            return c.get_branch_length() + (c.get_branch_length() + bl_augment * t);  });
    }

    double get_likelihood_for_diff_lambdas(const gene_family & gf, const clade * p_tree, const clade * p_lambda_tree, 
        int lambda_index, 
        std::vector<lambda*> & lambda_cache, 
        optimizer *opt,
        int max_root_family_size,
        int max_family_size)
    {
        const double bl_augment = 0.5;
        unique_ptr<clade> adjusted_tree(update_branchlength(p_tree, bl_augment, lambda_index));
        if (lambda_cache[lambda_index] == nullptr)
        {
            auto result = opt->optimize(optimizer_parameters());
            if (p_lambda_tree)
                lambda_cache[lambda_index] = new multiple_lambda(map<string, int>(), result.values);
            else
                lambda_cache[lambda_index] = new single_lambda(result.values[0]);
        }

        matrix_cache m(max_family_size + 1);
        m.precalculate_matrices(get_lambda_values(lambda_cache[lambda_index]), adjusted_tree->get_branch_lengths());
        auto probs = inference_prune(gf, m, lambda_cache[lambda_index], nullptr, adjusted_tree.get(), 1.0, max_root_family_size, max_family_size);
        return *max_element(probs.begin(), probs.end());
    }

    void compute_for_diff_lambdas_i(const user_data & data,
        std::vector<int> & lambda_index,
        std::vector<double> & pvalues,
        std::vector<lambda*> & lambda_cache,
        optimizer* p_opt
    )
    {
        auto references = build_reference_list(data.gene_families);

        matrix_cache cache(max(data.max_root_family_size, data.max_family_size) + 1);
        for (size_t i = 0; i < data.gene_families.size(); i += 1)
        {
            auto& pitem = data.gene_families[i];
            if (references[i] != i) continue;

            cache.precalculate_matrices(get_lambda_values(data.p_lambda), data.p_tree->get_branch_lengths());
            auto values = inference_prune(pitem, cache, data.p_lambda, data.p_error_model, data.p_tree, 1.0, data.max_root_family_size, data.max_family_size);
            double maxlh1 = *max_element(values.begin(), values.end());
            double prev = -1;
            double next = get_likelihood_for_diff_lambdas(pitem, data.p_tree, data.p_lambda_tree, 0, lambda_cache, p_opt, data.max_root_family_size, data.max_family_size);
            int j = 1;
            for (; prev < next; j++)
            {
                prev = next;
                next = get_likelihood_for_diff_lambdas(pitem, data.p_tree, data.p_lambda_tree, j, lambda_cache, p_opt, data.max_root_family_size, data.max_family_size);
            }
            pvalues[i] = (prev == maxlh1) ? 1 : 2 * (log(prev) - log(maxlh1));
            lambda_index[i] = j - 2;
        }
    }

    void likelihood_ratio_report(std::ostream & ost, const std::vector<gene_family> & families,
        const clade * pcafe,
        const std::vector<double> & pvalues,
        const std::vector<int> & plambda,
        const std::vector<lambda*> & lambda_cache)
    {
        for (size_t i = 0; i < families.size(); ++i)
        {
            ost << families[i].id() << "\t";
            pcafe->write_newick(ost, [](const clade* c) { return c->get_taxon_name(); });
            auto l = lambda_cache[plambda[i]];
            cout << "(" << plambda[i] << ", " << *l << ")\t" << pvalues[i] << "\t" << (pvalues[1] == 1 ? 1 : 1 - chi2cdf(pvalues[i], 1)) << endl;
        }
    }

    void lhr_for_diff_lambdas(const user_data & data, model *p_model)
    {
        std::vector<lambda*> lambda_cache(100);

        cout << "Running Likelihood Ratio Test 2....\n";

        std::vector<double> pvalues(data.gene_families.size());
        std::vector<int> lambdas(data.gene_families.size());

        auto lengths = data.p_tree->get_branch_lengths();
        auto longest_branch = *max_element(lengths.begin(), lengths.end());

        auto scorer = new lambda_optimizer(data.p_lambda, p_model, data.p_prior.get(), longest_branch, data.rootdist);

        optimizer opt(scorer);
        opt.quiet = true;
        compute_for_diff_lambdas_i(data, lambdas, pvalues, lambda_cache, &opt);

        likelihood_ratio_report(cout, data.gene_families, data.p_tree, pvalues, lambdas, lambda_cache);
    }
}