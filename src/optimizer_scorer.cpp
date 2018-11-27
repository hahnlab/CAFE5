#include <iomanip>
#include <iostream>

#include "optimizer_scorer.h"
#include "clade.h"
#include "lambda.h"
#include "base_model.h"
#include "gamma_core.h"
#include "gamma.h"

using namespace std;

double optimizer_scorer::calculate_score(double *values)
{
    prepare_calculation(values);

    if (!quiet)
    {
        report_precalculation();
    }

    _p_model->start_inference_processes(_p_lambda);

    double score = _p_model->infer_processes(_p_distribution);

    if (!quiet)
    {
        std::cout << "Score (-lnL): " << setw(15) << setprecision(14) << score << std::endl;
    }

    return score;
}


std::vector<double> lambda_optimizer::initial_guesses()
{
    branch_length_finder finder;
    _p_tree->apply_prefix_order(finder);
    std::vector<double> result(_p_lambda->count());
    for (auto& i : result)
    {
        i = 1.0 / finder.longest() * unifrnd();
    }
    return result;
}

void lambda_optimizer::prepare_calculation(double *values)
{
    _p_lambda->update(values);
}

void lambda_optimizer::report_precalculation()
{
    std::cout << "Lambda: " << *_p_lambda << std::endl;
}

std::vector<double> lambda_epsilon_optimizer::initial_guesses()
{
    std::vector<double> result(_p_lambda->count());
    for (auto& i : result)
    {
        i = 1.0 / _longest_branch * unifrnd();
    }

    current_guesses = _p_error_model->get_epsilons();
    result.insert(result.end(), current_guesses.begin(), current_guesses.end());

    return result;

}

void lambda_epsilon_optimizer::prepare_calculation(double *values)
{
    double * lambdas = values;
    double * epsilons = values + _p_lambda->count();

    _p_lambda->update(lambdas);
    map<double, double> replacements;
    for (size_t i = 0; i < current_guesses.size(); ++i)
    {
        replacements[current_guesses[i]] = epsilons[i];
        current_guesses[i] = epsilons[i];
    }

    _p_error_model->replace_epsilons(&replacements);
}

void lambda_epsilon_optimizer::report_precalculation()
{
    std::cout << "Calculating probability: epsilon=" << _p_error_model->get_epsilons().back()*2.0 << ", " << "lambda=" << *_p_lambda << std::endl;
}

void lambda_epsilon_optimizer::finalize(double *results)
{
    _p_lambda->update(results);
    _p_error_model->update_single_epsilon(results[_p_lambda->count()]);
}

gamma_optimizer::gamma_optimizer(gamma_model* p_model, root_equilibrium_distribution* prior) :
    optimizer_scorer(p_model->get_lambda(), p_model, prior),
    _p_gamma_model(p_model)
{

}

std::vector<double> gamma_optimizer::initial_guesses()
{
    return std::vector<double>({ unifrnd() });
}

void gamma_optimizer::prepare_calculation(double * values)
{
    double alpha = *values;
    _p_gamma_model->set_alpha(alpha);
}

void gamma_optimizer::report_precalculation()
{
    std::cout << "Attempting alpha: " << _p_gamma_model->get_alpha() << std::endl;
}

void gamma_optimizer::finalize(double * result)
{
    _p_gamma_model->set_alpha(*result);
}

void lambda_optimizer::finalize(double *results)
{
    _p_lambda->update(results);
}

gamma_lambda_optimizer::gamma_lambda_optimizer(const clade *p_tree, lambda *p_lambda, gamma_model * p_model, root_equilibrium_distribution *p_distribution) :
    optimizer_scorer(p_lambda, p_model, p_distribution),
    _p_tree(p_tree),
    _p_gamma_model(p_model)
{

}

std::vector<double> gamma_lambda_optimizer::initial_guesses()
{
    double alpha = unifrnd();

    std::vector<double> x(_p_gamma_model->get_gamma_cat_probs_count());
    std::vector<double> y(_p_gamma_model->get_lambda_multiplier_count());
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

void gamma_lambda_optimizer::prepare_calculation(double *values)
{
    _p_lambda->update(values);

    double alpha = values[_p_lambda->count()];
    _p_gamma_model->set_alpha(alpha);

    if (!quiet)
        _p_gamma_model->write_probabilities(cout);
}

void gamma_lambda_optimizer::report_precalculation()
{
    cout << "Attempting lambda: " << *_p_lambda << ", alpha: " << _p_gamma_model->get_alpha() << std::endl;
}

/// results consists of the desired number of lambdas and one alpha value
void gamma_lambda_optimizer::finalize(double *results) {
    _p_lambda->update(results);
    double alpha = results[_p_lambda->count()];
    _p_gamma_model->set_alpha(alpha);
}

