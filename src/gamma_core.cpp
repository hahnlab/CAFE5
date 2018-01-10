#include <assert.h>
#include <numeric>
#include <iomanip>
#include <cmath>

#include "gamma_core.h"
#include "gamma.h"
#include "process.h"
#include "root_equilibrium_distribution.h"
#include "reconstruction_process.h"


gamma_model::gamma_model(lambda* p_lambda, clade *p_tree, std::vector<gene_family>* p_gene_families, int max_family_size,
    int max_root_family_size, int n_gamma_cats, double fixed_alpha, std::map<int, int> *p_rootdist_map) :
    model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size) {
    if (p_rootdist_map != NULL)
        _rootdist_vec = vectorize_map(p_rootdist_map); // in vector form
    
    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
    _gamma_cats.resize(_rootdist_vec.size());
    set_alpha(fixed_alpha, _rootdist_vec.size());

    _total_n_families_sim = _rootdist_vec.size();
}

gamma_model::~gamma_model()
{
    for (size_t i = 0; i < _family_bundles.size(); ++i)
        _family_bundles[i].clear();
    _family_bundles.clear();
}

void gamma_model::print_results(std::ostream& ost)
{
    ost << "#FamilyID\tGamma Cat Median\tLikelihood of Category\tLikelihood of Family\tPosterior Probability\tSignificant" << endl;

    std::ostream_iterator<family_info_stash> out_it(ost, "\n");
    std::copy(results.begin(), results.end(), out_it);
}

//! Set alpha for gamma distribution
void gamma_model::set_alpha(double alpha, int n_families) {

    _alpha = alpha;
    if (_gamma_cat_probs.size() > 1)
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference

    vector<int>* cats = weighted_cat_draw(n_families, _gamma_cat_probs);
    _gamma_cats = *cats;
    delete cats;
}

void gamma_model::write_probabilities(ostream& ost)
{
    ost << "Gamma cat probs are: ";
    for (double d : _gamma_cat_probs)
        ost << d << ",";
    ost << endl;

    ost << "Lambda multipliers are: ";
    for (double d : _lambda_multipliers)
        ost << d << ",";
    ost << endl;
}

void gamma_model::start_inference_processes() {

    _family_bundles.clear();
    inference_process_factory factory(_ost, _p_lambda, _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec);
    for (auto i = _p_gene_families->begin(); i != _p_gene_families->end(); ++i)
    {
        factory.set_gene_family(&(*i));
        _family_bundles.push_back(gamma_bundle(factory, _lambda_multipliers));
    }
}

//! Populate _processes (vector of processes)
simulation_process* gamma_model::create_simulation_process(int family_number) {
    double lambda_bin = _gamma_cats[family_number];

    single_lambda *sl = dynamic_cast<single_lambda *>(_p_lambda);
    if (sl)
    {
        branch_length_finder lengths;
        _p_tree->apply_prefix_order(lengths);
        if (lengths.longest() * _lambda_multipliers[lambda_bin] * sl->get_single_lambda() > 1.0)
        {
            cerr << "WARNING: Probable saturation (branch length " << lengths.longest();
            cerr << " lambda " << sl->get_single_lambda() << " multiplier " << _lambda_multipliers[lambda_bin] << endl;
        }
    }

    return new simulation_process(_ost, _p_lambda, _lambda_multipliers[lambda_bin], _p_tree, _max_family_size, _max_root_family_size, _rootdist_vec, family_number); // if a single _lambda_multiplier, how do we do it?
}

std::vector<double> gamma_model::get_posterior_probabilities(std::vector<double> cat_likelihoods)
{
    size_t process_count = cat_likelihoods.size();

    vector<double> numerators(process_count);
    transform(cat_likelihoods.begin(), cat_likelihoods.end(), _gamma_cat_probs.begin(), numerators.begin(), multiplies<double>());

    double denominator = accumulate(numerators.begin(), numerators.end(), 0.0);
    vector<double> posterior_probabilities(process_count);
    transform(numerators.begin(), numerators.end(), posterior_probabilities.begin(), bind2nd(divides<double>(), denominator));

    return posterior_probabilities;
}

//! Infer bundle
double gamma_model::infer_processes(probability_calculator& calc, root_equilibrium_distribution *prior) {

    if (!_p_lambda->is_valid())
    {
        std::cout << "-lnL: " << log(0) << std::endl;
        return -log(0);
    }
    if (_alpha < 0)
    {
        std::cout << "-lnL: " << log(0) << std::endl;
        return -log(0);
    }

    using namespace std;
    initialize_rootdist_if_necessary();

    prior->initialize(_rootdist_vec);
    vector<double> all_bundles_likelihood(_family_bundles.size());

    bool success = true;
    calc.thread_cache();
    branch_length_finder lengths;
    _p_tree->apply_prefix_order(lengths);
    //_lambda_multipliers
    auto sl = dynamic_cast<single_lambda *>(_p_lambda);
    if (sl)
    {
        double lambda = sl->get_single_lambda();
        for (auto branch_length : lengths.result())
        {
            for (auto multiplier : _lambda_multipliers)
            {
                calc.get_matrix(_max_family_size+1, branch_length, lambda * multiplier);
            }
        }
    }

    vector<vector<family_info_stash>> pruning_results(_family_bundles.size());
#pragma omp parallel for
    for (int i = 0; i < _family_bundles.size(); ++i) {
        gamma_bundle& bundle = _family_bundles[i];

        try
        {
            vector<double> cat_likelihoods = bundle.prune(_gamma_cat_probs, prior, calc);

            double family_likelihood = accumulate(cat_likelihoods.begin(), cat_likelihoods.end(), 0.0);

            vector<double> posterior_probabilities = get_posterior_probabilities(cat_likelihoods);

            pruning_results[i].resize(cat_likelihoods.size());
            for (size_t k = 0; k < cat_likelihoods.size(); ++k)
            {
                pruning_results[i][k] = family_info_stash(i, bundle.get_lambda_likelihood(k), cat_likelihoods[k],
                    family_likelihood, posterior_probabilities[k], posterior_probabilities[k] > 0.95);
                //            cout << "Bundle " << i << " Process " << k << " family likelihood = " << family_likelihood << endl;
            }
            all_bundles_likelihood[i] = std::log(family_likelihood);
        }
        catch (runtime_error& ex)
        {
            // we got here because one of the gamma categories was saturated - reject this 
#pragma omp_critical
            success = false;
        }
    }
    calc.unthread_cache();
    if (!success)
        return -log(0);

    for (auto& stashes : pruning_results)
    {
        for (auto& stash : stashes)
        {
            results.push_back(stash);
        }
    }
    double final_likelihood = -accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 0.0);

    std::cout << "-lnL: " << final_likelihood << std::endl;

    return final_likelihood;
}

std::vector<double> gamma_model::initial_guesses()
{
    double alpha = unifrnd();

    std::vector<double> x(_gamma_cat_probs.size());
    std::vector<double> y(_lambda_multipliers.size());
    get_gamma(x, y, alpha); // passing vectors by reference

    double largest_multiplier = *max_element(y.begin(), y.end());
    branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);
    //double result = 1.0 / finder.result() * unifrnd();
    std::vector<double> lambdas(_p_lambda->count());
    const double longest_branch = finder.longest();
    generate(lambdas.begin(), lambdas.end(), [longest_branch, largest_multiplier] { return 1.0 / (longest_branch*largest_multiplier) * unifrnd(); });

    lambdas.push_back(alpha);
    return lambdas;
}

void gamma_model::set_current_guesses(double *guesses)
{
    _p_lambda->update(guesses);

    double alpha = guesses[_p_lambda->count()];
    set_alpha(alpha, _p_gene_families->size());

    cout << "Attempting lambda: " << *_p_lambda << ", alpha: " << alpha << std::endl;
}

void gamma_model::reconstruct_ancestral_states(probability_calculator *calc, root_equilibrium_distribution*prior)
{
    cout << "Reconstructing ancestral states using lambda = " << *_p_lambda << ", alpha = " << _alpha << endl;
    for (auto& bundle : _family_bundles)
    {
        bundle.set_values(calc, prior);
        bundle.reconstruct(_gamma_cat_probs);
    }
}

void gamma_model::print_reconstructed_states(std::ostream& ost)
{
    std::function<std::string()> f = [&]() {f = []() { return "\t"; }; return ""; };

    ost << "&Lambda multipliers: ";
    for (auto& i : _lambda_multipliers)
    {
        ost << f() << i;
    }
    ost << endl;

    auto rec = _family_bundles[0];
    auto order = rec.get_taxa();
    for (auto& it : order) {
        ost << "#" << it->get_taxon_name() << "\n";
    }

    for (auto& bundle : _family_bundles)
    {
        bundle.print_reconstruction(ost, order);
    }
}


