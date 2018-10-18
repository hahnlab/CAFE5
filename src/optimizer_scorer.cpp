#include <iomanip>
#include <iostream>

#include "optimizer_scorer.h"
#include "clade.h"
#include "lambda.h"
#include "base_model.h"
#include "gamma_core.h"

using namespace std;

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

double lambda_optimizer::calculate_score(double *values)
{
    _p_lambda->update(values);

    if (!quiet)
        std::cout << "Lambda: " << *_p_lambda << std::endl;

    _p_model->start_inference_processes();

    double score = _p_model->infer_processes(_p_distribution);

    if (!quiet)
        std::cout << "Score (-lnL): " << setw(15) << setprecision(14) << score << std::endl;

    return score;
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

double lambda_epsilon_optimizer::calculate_score(double *values)
{
    double * lambdas = values;
    double * epsilons = values + _p_lambda->count();

    if (*epsilons < 0 || *epsilons > .5)
    {
        return std::numeric_limits<double>::max();
    }

    _p_lambda->update(lambdas);
    map<double, double> replacements;
    for (size_t i = 0; i < current_guesses.size(); ++i)
    {
        replacements[current_guesses[i]] = epsilons[i];
        current_guesses[i] = epsilons[i];
    }

    _p_error_model->replace_epsilons(&replacements);

    if (!quiet)
    {
        std::cout << "Calculating probability: epsilon=" << _p_error_model->get_epsilons().back()*2.0 << ", " << "lambda=" << *_p_lambda << std::endl;
    }

    _p_model->start_inference_processes();

    double score = _p_model->infer_processes(_p_distribution);

    if (!quiet)
        std::cout << " Score with above error models: " << score << endl;

    return score;
}

void lambda_epsilon_optimizer::finalize(double *results)
{
    _p_lambda->update(results);
    _p_error_model->update_single_epsilon(results[_p_lambda->count()]);
}

/// results consists of the desired number of lambdas and one alpha value
void gamma_lambda_optimizer::finalize(double *results) {
    _p_lambda->update(results);
    double alpha = results[_p_lambda->count()];
    _p_model->set_alpha(alpha, _p_model->get_gene_family_count());
}

std::vector<double> gamma_optimizer::initial_guesses()
{
    return std::vector<double>({ unifrnd() });
}

double gamma_optimizer::calculate_score(double * values)
{
    double alpha = *values;
    _p_model->set_alpha(alpha, _p_model->get_gene_family_count());

    std::cout << "Attempting alpha: " << alpha << std::endl;

    _p_model->start_inference_processes();

    return _p_model->infer_processes(_p_distribution);
}

void gamma_optimizer::finalize(double * result)
{
    _p_model->set_alpha(*result, _p_model->get_gene_family_count());
}

void lambda_optimizer::finalize(double *results)
{
    _p_lambda->update(results);
}
