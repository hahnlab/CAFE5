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
    
    this->initialize_with_alpha(n_gamma_cats, _rootdist_vec.size(), fixed_alpha);
    _total_n_families_sim = _rootdist_vec.size();
}

//! Simulation: gamma_model constructor when just alpha is provided.
gamma_model::gamma_model(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
    int n_gamma_cats, double alpha) : model(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
    _gamma_cat_probs(n_gamma_cats), _lambda_multipliers(n_gamma_cats) {

    if (!rootdist_vec.empty()) {
        _rootdist_bins.push_back(rootdist_vec); // just 1st element
    }

    else {
        _rootdist_vec = uniform_dist(total_n_families, 1, max_family_size); // the user did not specify one... using uniform from 1 to max_family_size!
        _rootdist_bins.push_back(_rootdist_vec); // just 1st element (but we could specify different root dists for each lambda bin)
    }

    if (n_gamma_cats > 1) {
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference
        auto cats = weighted_cat_draw(total_n_families, _gamma_cat_probs);
        _gamma_cats = *cats;
        delete cats;
    }

    for (auto i = _gamma_cat_probs.begin(); i != _gamma_cat_probs.end(); ++i) {
        cout << "Should be all the same probability: " << *i << endl;
    }

    for (auto i = _lambda_multipliers.begin(); i != _lambda_multipliers.end(); ++i) {
        cout << "Lambda multiplier (rho) is: " << *i << endl;
    }

    for (auto i = _gamma_cats.begin(); i != _gamma_cats.end(); ++i) {
        cout << "Gamma category is: " << *i << endl;
    }

}

//! Simulation: model constructor when alpha is not provided and membership is provided directly, in addition to the lambda multipliers (something the user decides to use)
gamma_model::gamma_model(ostream & ost, lambda* lambda, clade *p_tree, int max_family_size, int total_n_families, vector<int> rootdist_vec,
    vector<int>& cats, vector<double>&mul) : model(ost, lambda, p_tree, max_family_size, total_n_families, rootdist_vec),
    _gamma_cats(cats), _lambda_multipliers(mul) {}

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


void gamma_model::initialize_with_alpha(int n_gamma_cats, int n_families, double alpha)
{
    adjust_n_gamma_cats(n_gamma_cats);
    adjust_family_gamma_membership(n_families);
    set_alpha(alpha, n_families);
}

void gamma_model::initialize_without_alpha(int n_gamma_cats, int n_families, vector<double> lambda_multipliers, std::vector<int> gamma_cats)
{
    assert(lambda_multipliers.size() == n_gamma_cats);
    assert(gamma_cats.size() == n_families);
    //    adjust_n_gamma_cats(n_gamma_cats);
    //    adjust_family_gamma_membership(n_families);
    //    set_alpha(alpha, n_families);

    set_lambda_multipliers(lambda_multipliers);
    _gamma_cats = gamma_cats;
}

//! Resize all gamma-related vectors according to provided number (integer) of gamma categories
void gamma_model::adjust_n_gamma_cats(int n_gamma_cats) {
    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
}

//! Resize gamma_cats vector that assigns gamma class membership of families to be inferred/simulated
void gamma_model::adjust_family_gamma_membership(int n_families) {
    _gamma_cats.resize(n_families);
}

//! Set alpha for gamma distribution
void gamma_model::set_alpha(double alpha, int n_families) {
    _alpha = alpha;
    if (_gamma_cats.size() > 1)
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference

    vector<int>* cats = weighted_cat_draw(n_families, _gamma_cat_probs);
    _gamma_cats = *cats;
    delete cats;

//    for (std::vector<double>::iterator it = _gamma_cat_probs.begin(); it != _gamma_cat_probs.end(); ++it) {
//        cout << "Gamma cat prob is : " << *it << endl;
//    }
}

//! Set lambda multipliers for each gamma category
void gamma_model::set_lambda_multipliers(std::vector<double> lambda_multipliers) {
    _lambda_multipliers = lambda_multipliers;
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
double gamma_model::infer_processes(root_equilibrium_distribution *prior) {

    using namespace std;
    initialize_rootdist_if_necessary();

    prior->initialize(_rootdist_vec);
    vector<double> all_bundles_likelihood(_family_bundles.size());

    for (int i = 0; i < _family_bundles.size(); ++i) {
//        cout << endl << "About to prune a gamma bundle." << endl;
        gamma_bundle& bundle = _family_bundles[i];

        vector<double> cat_likelihoods = bundle.prune(_gamma_cat_probs, prior);
        double family_likelihood = accumulate(cat_likelihoods.begin(), cat_likelihoods.end(), 0.0);

        vector<double> posterior_probabilities = get_posterior_probabilities(cat_likelihoods);

        for (size_t k = 0; k < cat_likelihoods.size(); ++k)
        {            
            results.push_back(family_info_stash(i, bundle.get_lambda_likelihood(k), cat_likelihoods[k], 
                family_likelihood, posterior_probabilities[k], posterior_probabilities[k] > 0.95));
//            cout << "Bundle " << i << " Process " << k << " family likelihood = " << family_likelihood << endl;
        }

        all_bundles_likelihood[i] = std::log(family_likelihood);
//        cout << "Bundle " << i << " family likelihood = " << family_likelihood << endl;

//        cout << "Likelihood of family " << i << " = " << all_bundles_likelihood[i] << endl;
    }

    double final_likelihood = -accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 0.0);

    std::cout << "-lnL: " << final_likelihood << std::endl;

    return final_likelihood;
}

std::vector<double> gamma_model::initial_guesses()
{
    max_branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);
    double result = 1.0 / finder.result() * unifrnd();
    std::vector<double> lambdas(_p_lambda->count());
    for (auto& i : lambdas)
    {
        i = 1.0 / finder.result() * unifrnd();
    }

    cout << "Initial lambda: " << result << std::endl;
    double alpha = unifrnd();
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
    for (auto& bundle : _family_bundles)
    {
        bundle.set_values(calc, prior);
        bundle.reconstruct(_gamma_cat_probs);
    }
}

void gamma_model::print_reconstructed_states(std::ostream& ost)
{
    cerr << "gamma_model::print_reconstructed_states\n";

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


