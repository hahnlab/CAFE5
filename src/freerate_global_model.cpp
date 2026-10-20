#include <algorithm>
#include <numeric>
#include <iomanip>
#include <random>

#include <boost/math/tools/minima.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "freerate_global_model.h"
#include "user_data.h"
#include "inference_pruner.h"
#include "matrix_cache.h"
#include "prior.h"
#include "lambda.h"
#include "optimizer_scorer.h"   // for definition of compute_distribution_mean
#include "spearman.h"
#include "fitch_margoliash.h"

using namespace std;
using namespace Eigen;

extern std::mt19937 randomizer_engine; // seeding random number engine
typedef vector<optional_probabilities> ASV;

const clade * get_sibling(const clade *c)
{
    if (c->is_root())
        return nullptr;
    auto p = c->get_parent();
    auto sib = find_if(p->descendant_begin(), p->descendant_end(), [&](clade *d) {
        return d != c;
    });
    return sib == p->descendant_end() ? nullptr : *sib;
}

double get_log_likelihood(const ASV& partial_likelihoods, const std::vector<double>& priors)
{
    std::vector<double> all_transcripts_likelihood(partial_likelihoods.size());

    transform(partial_likelihoods.begin(), partial_likelihoods.end(), all_transcripts_likelihood.begin(), [&](const optional_probabilities& a) {
        vector<double> t(a.probabilities().size());
        copy(a.probabilities().begin(), a.probabilities().end(), t.begin());
        return log(compute_prior_likelihood(t, priors));
    });

    double result = std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);
    if (std::isinf(result))
        result = 10000;
    return result;
}

void init_tips(const user_data &ud, transcript_clade_map<optional_probabilities> &probs, const matrix_cache &cache)
{
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade *p)
             {
        if (p->is_leaf())
        {
            for (auto& transcript : ud.gene_transcripts)
            {
                probs.get(&transcript, p).reserve(cache.create_vector());
                probs.get(&transcript, p).initialize(transcript.get_expression_value(p->get_taxon_name()), ud.bounds);
            }
        } });
}


freerate_global_model::freerate_global_model(bool values_are_ratios, optimizer_parameters params) :
    model(nullptr, nullptr, nullptr), _values_are_ratios(values_are_ratios), _optimizer_params(params)
{
}

// This has a lot of elements in common with compute_node_probability. Can they be refactored?
double compute_node_likelihood(const clade* d,
                    transcript_clade_map<optional_probabilities>& probs, clademap<std::pair<double, double>>& sigmas,
                    const vector<double>& priors_by_bin, double sigsqd, matrix_cache& cache, const user_data& ud)
{
    auto p = d->get_parent();
    auto sib = get_sibling(d);
    cache.precalculate_matrices({ sigsqd, sigmas[sib].first }, { d->get_branch_length(), sib->get_branch_length() }, ud.bounds);

    for (auto& transcript : ud.gene_transcripts)
    {
        probs.get(&transcript, p).reserve(cache.create_vector());
    }

#pragma omp parallel for
    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i)
    {
        auto& transcript = ud.gene_transcripts[i];
        auto& node_probs = probs.get(&transcript, d);
        auto& parent_probs = probs.get(&transcript, p);
        auto& sib_probs = probs.get(&transcript, sib);
        if (node_probs.hasValue() || sib_probs.hasValue())
        {
            parent_probs.setOne();
        }

        if (node_probs.hasValue())
        {
            const MatrixXd& m = cache.get_matrix(d->get_branch_length(), sigsqd);
            VectorXd result = m * node_probs.probabilities();
            parent_probs.multiply_elements(result);
        }
        if (sib_probs.hasValue())
        {
            const MatrixXd& m = cache.get_matrix(sib->get_branch_length(), sigmas[sib].first);
            VectorXd result = m * sib_probs.probabilities();
            parent_probs.multiply_elements(result);
        }
    }

    auto lnl = get_log_likelihood(probs.get_all(p), priors_by_bin);
    LOG(DEBUG) << "Node " << d->get_ape_index() << " Sigma^2:" << sigsqd << " Score (-lnL): " << lnl;
    return -lnl;
}

/**
 * @brief Returns the maximum value to be used for sigma^2 optimization.
 *
 * This function iterates through the provided map of sigma values for each clade,
 * finds the maximum sigma^2 value, and returns either 100 times this value or 20.0,
 * whichever is smaller. This is used to set an upper bound for the optimizer when
 * searching for the optimal sigma^2 for a branch.
 *
 * @param sigmas A map from clade pointers to pairs of (sigma^2, auxiliary value).
 * @return The maximum allowed value for sigma^2 during optimization.
 */
double get_optimizer_max(clademap<std::pair<double, double>>& sigmas)
{
    const double default_max = 20.0;
    double max_sigma = 0.0;
    for (const auto& kv : sigmas) {
        if (kv.second.first > max_sigma) {
            max_sigma = kv.second.first;
        }
    }
    if (max_sigma == 0.0) {
        return default_max;
    }
    return std::min(max_sigma * 100.0, default_max);
}

/* Init all branches with a sigma
To calculate AF:
        ASV[B] = get_ASV(childB branch length, sigma[B], childB init)
        optimize ss:
                ASV[A] = get_ASV(childA branch length, ss, childA init)
                ASV[F] = ASV[AF] * ASV[BF]
                likelihood = ASV[F] * prior


To calculate BF:
        ASV[BF] = optimized value of get_ASV(childA bl, childB bl, sigma, ASV[AF], childB init)
To calculate GF:
        ASV[GC] = get_ASV(childC branch length, childF branch length, sigma_init, childC init, ASV[AF] * ASV[BF])
        ASV[GF] = optimized value of get_ASV(childC bl, childF bl, sigma, ASV[GC], childC init)
 */
double freerate_global_model::optimize_sigmas(const user_data& ud, const clademap<prior>& priors)
{
    matrix_cache cache;

    transcript_clade_map<optional_probabilities> probs;

    const double optimizer_max = get_optimizer_max(_sigmas);
    init_tips(ud, probs, cache);
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* p) {
        if (!p->is_leaf())
        {
            auto priors_by_bin = get_priors(cache, ud.bounds, &priors.at(p));
            cladevector descendants(p->descendant_begin(), p->descendant_end());
            shuffle(descendants.begin(), descendants.end(), randomizer_engine);  
            for (auto d : descendants)
            {
                using boost::math::tools::brent_find_minima;

                LOG(DEBUG) << "Optimizing node " << d->get_ape_index() << " with sibling sigma^2 " << _sigmas[get_sibling(d)].first;
                std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                    return compute_node_likelihood(d, probs, _sigmas, priors_by_bin, sigsqd, cache, ud);
                }, 0.0, optimizer_max, 15);
                _sigmas[d] = r;
                LOG(INFO) << "Node " << d->get_ape_index() << " Sigma^2:" << r.first << " Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;
            }
        }
    });

    auto root_priors = get_priors(cache, ud.bounds, &priors.at(ud.p_tree));
    auto result = -get_log_likelihood(probs.get_all(ud.p_tree), root_priors);
    LOG(INFO) << "Final tree score (-lnL): " << std::setw(15) << std::setprecision(14) << result;
    return result;
}

fitch_margoliash::distance_matrix invert_matrix(const fitch_margoliash::distance_matrix& distance_matrix)
{
    fitch_margoliash::distance_matrix inverted;
    for (const auto& kv : distance_matrix) {
        inverted[std::make_pair(kv.first.first, kv.first.second)] = 1.0 - kv.second;
    }
    return inverted;
}

void freerate_global_model::initialize_sigmas(const clade* p_tree, const transcript_vector& transcripts)
{
    unique_ptr<clade> p_weight_tree;
    if (transcripts.size() > 1)
    {
        auto correlation = spearman_correlation_by_species(transcripts);
        auto distance = invert_matrix(correlation);
        fitch_margoliash fm(distance, transcripts[0].get_species());

        auto p = fm.build_tree();
        p->normalize();
        p_weight_tree.reset(p);
    }

    double distmean = compute_distribution_mean(p_tree, transcripts);
    LOG(DEBUG) << "Distribution mean: " << distmean;
    
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&](const clade* p) {
        if (p->is_root() || !p_weight_tree)
            _sigmas[p] = pair<double,double>(distmean,-1);
        else 
        {
            auto weight = 0.0;
            auto matched_clade = p_weight_tree->find_descendant(p->get_taxon_name());
            if (matched_clade)
            {   
                weight = matched_clade->get_branch_length();
                LOG(DEBUG) << "Found a weight for " << p->get_taxon_name() << ". It is " << weight;
            }
            else
            {
                vector<string> descendants(distance(p->descendant_begin(), p->descendant_end()));
                transform(p->descendant_begin(), p->descendant_end(), descendants.begin(), [](const clade* c) { return c->get_taxon_name(); });
                auto d1 = p_weight_tree->find_descendant(descendants[0]);
                auto d2 = p_weight_tree->find_descendant(descendants[1]);
                weight = distance(d1, d2);

                LOG(DEBUG) << "There is no weight for " << p->get_taxon_name() << ". Computed distance between " << descendants[0] << " and " << descendants[1] << " is " << weight;
            }

            weight = distmean * weight;
            _sigmas[p] =  pair<double,double>(weight,-1);  
        }
        string name = p->is_leaf() ? " (" + p->get_taxon_name() + ")" : "";
        LOG(DEBUG) << "Initial sigma for Node " << p->get_ape_index() << name << " set to " << _sigmas[p].first;

    });
    LOG(DEBUG) << "Initial sigmas set to " << distmean;
}

double freerate_global_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma)
{
    _root_ape_index = ud.p_tree->get_ape_index();

    auto priors = compute_tree_priors(ud.p_tree, ud.gene_transcripts, _values_are_ratios);
    for (auto& a : priors)
    {
        LOG(DEBUG) << "Node " << a.first->get_ape_index() << " Prior: " << a.second;
    }
    initialize_sigmas(ud.p_tree, ud.gene_transcripts);

    double score = 100;
    for (int i = 0; i<_optimizer_params.max_iterations; ++i)
    {
        LOG(INFO) << "Iteration " << i+1;
        _monitor.Event_InferenceAttempt_Started();
        double lnl = optimize_sigmas(ud, priors);
        if (abs(score - lnl) < 10)
            break;
        score = lnl;
    }
    return score;
}

sigma_optimizer_scorer* freerate_global_model::get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups)
{
    return nullptr;
}

class freerate_global_reconstruction : public reconstruction
{
public:
    freerate_global_reconstruction(const user_data& ud, matrix_cache* p_calc) : reconstruction(ud.gene_transcripts) {};
    virtual ~freerate_global_reconstruction() {};

    transcript_clade_map<node_reconstruction> _reconstructions;

    virtual node_reconstruction get_internal_node_value(const gene_transcript& transcript, const clade* c) const override
    {
        if (!_reconstructions.find(&transcript, c))
            throw missing_expression_value(transcript.id(), "");

        return _reconstructions.get(&transcript, c);
    }
};

void compute_transcript_likelihoods(const clade* d, const gene_transcript& transcript, const matrix_cache& cache, transcript_clade_map<optional_probabilities>& probs, const clademap<std::pair<double, double>>& sigmas)
{
    if (d->is_root())
        return;

    auto p = d->get_parent();
    probs.get(&transcript, p).reserve(cache.create_vector());

    auto sib = get_sibling(d);

    auto& node_probs = probs.get(&transcript, d);
    auto& parent_probs = probs.get(&transcript, p);
    auto& sib_probs = probs.get(&transcript, sib);

    if (node_probs.hasValue() || sib_probs.hasValue())
    {
        parent_probs.setOne();
    }

    if (node_probs.hasValue())
    {
        const MatrixXd& m = cache.get_matrix(d->get_branch_length(), sigmas.at(d).first);
        VectorXd result = m * node_probs.probabilities();
        parent_probs.multiply_elements(result);
    }
    if (sib_probs.hasValue())
    {
        const MatrixXd& m = cache.get_matrix(sib->get_branch_length(), sigmas.at(sib).first);
        VectorXd result = m * sib_probs.probabilities();
        parent_probs.multiply_elements(result);
    }
}

extern node_reconstruction get_value(const VectorXd& likelihood, boundaries bounds);

reconstruction* freerate_global_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    std::vector<double> sigma_values;
    sigma_values.reserve(_sigmas.size());
    for (const auto& kv : _sigmas) {
        sigma_values.push_back(kv.second.first);
    }
    p_calc->precalculate_matrices(sigma_values, ud.p_tree->get_branch_lengths(), ud.bounds);

    auto result = new freerate_global_reconstruction(ud, p_calc);

    transcript_clade_map<optional_probabilities> probs;

    init_tips(ud, probs, *p_calc);

    for (auto& transcript : ud.gene_transcripts)
    {
        for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* c) {
            compute_transcript_likelihoods(c, transcript, *p_calc, probs, _sigmas);
            auto& p = probs.get(&transcript, c);
            if (!c->is_leaf() && p.hasValue())
                result->_reconstructions.set(&transcript, c, get_value(p.probabilities(), ud.bounds));
        });
    }

    return result;
}

void freerate_global_model::write_extra_vital_statistics(std::ostream& ost)
{
    ost << "Computed sigma2 by node:\n";
    vector<pair<int,double>> sorted;
    transform(_sigmas.begin(), _sigmas.end(), back_inserter(sorted), [](const pair<const clade*, pair<double, double>>& a) { return make_pair(a.first->get_ape_index(), a.second.first); });
    sort(sorted.begin(), sorted.end(), [](const pair<int, double>& a, const pair<int, double>& b) { return a.first < b.first; });
    for (auto& a : sorted)
    {
        if (a.first == _root_ape_index)
            continue;

        ost << a.first << ":" << a.second << "\n";
    }

    if (!sorted.empty()) {
        // Find minimum value
        double min_sigma = sorted[0].second;
        for (const auto& a : sorted) {
            if (a.second < min_sigma) min_sigma = a.second;
        }
        // Count occurrences
        int count = 0;
        for (const auto& a : sorted) {
            if (a.second == min_sigma) ++count;
        }
        ost << "Minimum sigma2: " << min_sigma << " occurred " << count << " times\n";
    }
}



TEST_CASE("freerate_model")
{
    new freerate_global_model(false, optimizer_parameters());
}

TEST_CASE("freerate_model optimizes a branch length")
{
    freerate_global_model m(false, optimizer_parameters());

    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("uniform", 0.0, 1.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    m.infer_transcript_likelihoods(ud, nullptr);
}

TEST_CASE("get_log_likelihood")
{
    ASV asv(1);
    Eigen::VectorXd v(5);
    v.setOnes();
    asv[0].set(v);
    auto lnl = get_log_likelihood(asv, vector<double>({ 1.0, 1.0, 1.0, 1.0, 1.0 }));
    CHECK_EQ(doctest::Approx(1.60944), lnl);
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("write_extra_vital_statistics writes a sigma for each internal node")
{
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("gamma", 0.375, 1600.0);
    ud.p_tree = parse_newick("((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;");

    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 0.5);
    ud.gene_transcripts[0].set_expression_value("B", 1);
    ud.gene_transcripts[0].set_expression_value("C", 1.5);
    ud.gene_transcripts[0].set_expression_value("D", 2);
    ud.gene_transcripts[0].set_expression_value("E", 2.5);

    freerate_global_model m(false, optimizer_parameters());
    m.initialize_sigmas(ud.p_tree, ud.gene_transcripts);

    ostringstream ost;
    m.write_extra_vital_statistics(ost);
    CHECK_STREAM_CONTAINS(ost, "Computed sigma2 by node:");
    CHECK_STREAM_CONTAINS(ost, "7:0.428082");
    CHECK_STREAM_CONTAINS(ost, "8:0.428082");
    CHECK_STREAM_CONTAINS(ost, "9:0.428082");
}

TEST_CASE("get_parent_likelihood")
{
    user_data ud;
    ud.p_sigma = nullptr;
    ud.p_prior = new prior("uniform", 0.0, 100.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    matrix_cache cache;
    transcript_clade_map<optional_probabilities> probs;
    const clade* p = ud.p_tree;
    const clade* sib = get_sibling(*p->descendant_begin());
    const clade* d = *p->descendant_begin();
    auto transcript = &ud.gene_transcripts[0];
    probs.get(transcript, d).reserve(cache.create_vector());
    probs.get(transcript, d).initialize(transcript->get_expression_value("A"), ud.bounds);
    probs.get(transcript, sib).reserve(cache.create_vector());
    probs.get(transcript, sib).initialize(transcript->get_expression_value("B"), ud.bounds);

    vector<double> priors_by_bin = get_priors(cache, ud.bounds, ud.p_prior);
    double sigsqd = 1.0;
    clademap<std::pair<double, double>> sigmas;
    sigmas[sib]= pair<double,double>(1.5, -1.0);

    double likelihood = compute_node_likelihood(d, probs, sigmas, priors_by_bin, sigsqd, cache, ud);

    CHECK_EQ(doctest::Approx(8.07594), likelihood);
}

TEST_CASE("compute_transcript_likelihoods")
{
    // Setup user_data and tree
    user_data ud;
    ud.p_sigma = nullptr;
    ud.p_prior = new prior("uniform", 0.0, 10.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);

    // Add a gene transcript with values for A and B
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1.0);
    ud.gene_transcripts[0].set_expression_value("B", 2.0);

    // Prepare matrix_cache and ASV probabilities
    matrix_cache cache;
    cache.precalculate_matrices({1,1}, ud.p_tree->get_branch_lengths(), ud.bounds);
    transcript_clade_map<optional_probabilities> probs;
    const clade* p = ud.p_tree;
    const clade* d = *p->descendant_begin();
    const clade* sib = get_sibling(d);

    for (const auto& transcript : ud.gene_transcripts)
    {
        probs.get(&transcript, d).reserve(cache.create_vector());
        probs.get(&transcript, d).initialize(transcript.get_expression_value("A"), ud.bounds);
        probs.get(&transcript, sib).reserve(cache.create_vector());
        probs.get(&transcript, sib).initialize(transcript.get_expression_value("B"), ud.bounds);
    }
    
    // Setup sigmas for both children
    clademap<std::pair<double, double>> sigmas;
    sigmas[d] = std::make_pair(1.0, -1.0);
    sigmas[sib] = std::make_pair(1.0, -1.0);

    // Compute transcript likelihoods for gene 0
    auto t = &ud.gene_transcripts[0];
    compute_transcript_likelihoods(d, *t, cache, probs, sigmas);

    // Check that parent probabilities are set (should not be empty)
    CHECK(probs.get(t,p).hasValue());
    // Optionally, check that the probabilities vector is not all zeros
    auto v = probs.get(t,p).probabilities();
    bool nonzero = false;
    for (int i = 0; i < v.size(); ++i) {
        if (v[i] != 0.0) {
            nonzero = true;
            break;
        }
    }
    CHECK(nonzero);

}

TEST_CASE("reconstruct_ancestral_states")
{
    user_data ud;
    ud.p_tree = parse_newick("((A:1,B:1):1,(C:1,D:1):1);");
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);
    ud.gene_transcripts[0].set_expression_value("C", 5);
    ud.gene_transcripts[0].set_expression_value("D", 6);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily2", "", ""));
    ud.gene_transcripts[1].set_expression_value("A", 3);
    ud.gene_transcripts[1].set_expression_value("B", 4);
    ud.gene_transcripts[1].set_expression_value("C", 7);
    ud.gene_transcripts[1].set_expression_value("D", 8);
    ud.bounds = boundaries(0, 20);
    sigma_squared sl(0.05);
    ud.p_sigma = &sl;

    freerate_global_model model(false, optimizer_parameters());
    model.initialize_sigmas(ud.p_tree, ud.gene_transcripts);

    matrix_cache calc;

    std::unique_ptr<freerate_global_reconstruction> rec(
        dynamic_cast<freerate_global_reconstruction*>(model.reconstruct_ancestral_states(ud, &calc))
    );

    // Two families * threee internal nodes = 6 reconstructions
    CHECK_EQ(6, rec->_reconstructions.count());
}

TEST_CASE("get_optimizer_max returns 20 if 100 * largest sigma is above that")
{
    clade a,b,c;
    clademap<std::pair<double, double>> sigmas;
    sigmas[&a] = std::make_pair(0.2, -1);
    sigmas[&b] = std::make_pair(0.3, -1);
    sigmas[&c] = std::make_pair(0.4, -1);

    CHECK_EQ(doctest::Approx(20.0), get_optimizer_max(sigmas));
}

TEST_CASE("get_optimizer_max returns 100 * largest sigma if that is below 20")
{
    clade a,b,c;
    clademap<std::pair<double, double>> sigmas;
    sigmas[&a] = std::make_pair(0.02, -1);
    sigmas[&b] = std::make_pair(0.03, -1);
    sigmas[&c] = std::make_pair(0.04, -1);

    CHECK_EQ(doctest::Approx(4.0), get_optimizer_max(sigmas));
}

TEST_CASE("invert_matrix")
{
    fitch_margoliash::distance_matrix matrix;
    matrix[std::make_pair("A", "B")] = 0.1;
    matrix[std::make_pair("A", "C")] = 0.2;
    matrix[std::make_pair("B", "C")] = 0.3;

    fitch_margoliash::distance_matrix inverted = invert_matrix(matrix);

    CHECK_EQ(inverted[std::make_pair("A", "B")], 0.9);
    CHECK_EQ(inverted[std::make_pair("A", "C")], 0.8);
    CHECK_EQ(inverted[std::make_pair("B", "C")], 0.7);
}

TEST_CASE("initialize_sigmas sets initial sigma values based on distance matrix")
{
    unique_ptr<clade> tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
    vector<gene_transcript> transcripts;

    gene_transcript t1("Gene1", "", "");
    t1.set_expression_value("A", 1.0);
    t1.set_expression_value("B", 2.0);
    t1.set_expression_value("C", 3.0);
    t1.set_expression_value("D", 4.0);
    transcripts.push_back(t1);  

    gene_transcript t2("Gene2", "", "");
    t2.set_expression_value("A", 2.0);
    t2.set_expression_value("B", 3.0);
    t2.set_expression_value("C", 0.5);
    t2.set_expression_value("D", 2.5);
    transcripts.push_back(t2);

    freerate_global_model m(false, optimizer_parameters());
    m.initialize_sigmas(tree.get(), transcripts);

    ostringstream ost;
    m.write_extra_vital_statistics(ost);
    CHECK_STREAM_CONTAINS(ost, "Computed sigma2 by node:");
    CHECK_STREAM_CONTAINS(ost, "5:0.708333");
    CHECK_STREAM_CONTAINS(ost, "6:4.95833");
    CHECK_STREAM_CONTAINS(ost, "7:4.95833");

}s