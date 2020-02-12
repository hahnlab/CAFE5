class clade;
class gene_family;
class user_data;
class lambda;
class optimizer_scorer;
class model;

#include <vector>

namespace LikelihoodRatioTest {
    clade* update_branchlength(const clade* lambda_tree, double bl_augment, int t);

    double get_likelihood_for_diff_lambdas(const gene_family& gf, const clade* p_tree,
        const clade* p_lambda_tree, int t, std::vector<lambda*>& lambda_cache, optimizer* p_opt, int, int);

    void compute_for_diff_lambdas_i(const user_data& data,
        std::vector<int>& lambda_index,
        std::vector<double>& pvalues,
        std::vector<lambda*>& lambda_cache,
        optimizer* p_opt
    );

    void lhr_for_diff_lambdas(const user_data& data, model *p_model);
}

