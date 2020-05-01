#include <iomanip>
#include <iostream>
#include <random>

#include "easylogging++.h"

#include "optimizer_scorer.h"
#include "clade.h"
#include "lambda.h"
#include "base_model.h"
#include "gamma_core.h"
#include "gamma.h"
#include "error_model.h"
#include "root_equilibrium_distribution.h"

#define GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA 1.75

extern std::mt19937 randomizer_engine;

using namespace std;

double inference_optimizer_scorer::calculate_score(const double *values)
{
    prepare_calculation(values);

    if (!quiet)
    {
        report_precalculation();
    }

    double score = _p_model->infer_family_likelihoods(*_p_distribution, _p_lambda);

    if (std::isnan(score)) score = -log(0);

    return score;
}

//Inititial Guess multiplies the 1/longest branch by a random draw from a normal 
//distribution centered such that it will start around a value for lambda of 0.002
std::vector<double> lambda_optimizer::initial_guesses()
{
    double distmean = 0.002/(1.0 / _longest_branch);
    std::vector<double> result(_p_lambda->count());
    //std::uniform_real_distribution<double> distribution(0.0, 1.0); Insert a prior distribution (above to start from a biologically realistic rate)
    std::normal_distribution<double> distribution(distmean,0.2);
    for (auto& i : result)
    {
    	i=1.0 / _longest_branch * distribution(randomizer_engine);
    	while (i<0)
    	{
    		i=1.0 / _longest_branch * distribution(randomizer_engine);
    	}
    }
    return result;
}

void lambda_optimizer::prepare_calculation(const double *values)
{
    _p_lambda->update(values);
}

void lambda_optimizer::report_precalculation()
{
    LOG(INFO) << "Lambda: " << *_p_lambda;
}

void lambda_optimizer::finalize(double *results)
{
    _p_lambda->update(results);
}

std::vector<double> lambda_epsilon_optimizer::initial_guesses()
{
    auto result = _lambda_optimizer.initial_guesses();

    current_guesses = _p_error_model->get_epsilons();
    result.insert(result.end(), current_guesses.begin(), current_guesses.end());

    return result;

}

void lambda_epsilon_optimizer::prepare_calculation(const double *values)
{
    auto lambdas = values;
    auto epsilons = values + _p_lambda->count();

    _lambda_optimizer.prepare_calculation(lambdas);

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
    LOG(INFO) << "Calculating probability: epsilon=" << _p_error_model->get_epsilons().back()*2.0 << ", " << "lambda=" << *_p_lambda;
}

void lambda_epsilon_optimizer::finalize(double *results)
{
    _lambda_optimizer.finalize(results);
    _p_error_model->update_single_epsilon(results[_p_lambda->count()]);
}

gamma_optimizer::gamma_optimizer(gamma_model* p_model, const root_equilibrium_distribution* prior) :
    inference_optimizer_scorer(p_model->get_lambda(), p_model, prior),
    _p_gamma_model(p_model)
{

}
//In the first line that is commented out, Alpha is initiated by randomly drawing from an exponential distribution with a mean of 1.75
//It seems a gamma distribution with an alpha of 4 and a beta scaling factor of 0.25 works better. It has a mean of 1 and 70% of the density is between .5 and 1.5
std::vector<double> gamma_optimizer::initial_guesses()
{
    //std::exponential_distribution<double> distribution(GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA);
    std::gamma_distribution<double> distribution(4.0,0.25);
    return std::vector<double>({ distribution(randomizer_engine) });
}

void gamma_optimizer::prepare_calculation(const double * values)
{
    double alpha = *values;
    _p_gamma_model->set_alpha(alpha);
}

void gamma_optimizer::report_precalculation()
{
    LOG(INFO) << "Attempting alpha: " << _p_gamma_model->get_alpha();
}

void gamma_optimizer::finalize(double * result)
{
    _p_gamma_model->set_alpha(*result);
}

double gamma_optimizer::get_alpha() const
{
    return _p_gamma_model->get_alpha();
}

gamma_lambda_optimizer::gamma_lambda_optimizer(lambda *p_lambda, gamma_model * p_model, const root_equilibrium_distribution *p_distribution, double longest_branch) :
    inference_optimizer_scorer(p_lambda, p_model, p_distribution),
    _lambda_optimizer(p_lambda, p_model, p_distribution, longest_branch),
    _gamma_optimizer(p_model, p_distribution)
{
}

std::vector<double> gamma_lambda_optimizer::initial_guesses()
{
    auto values = _lambda_optimizer.initial_guesses();
    auto alpha = _gamma_optimizer.initial_guesses();

    values.insert(values.end(), alpha.begin(), alpha.end());
    return values;

}

void gamma_lambda_optimizer::prepare_calculation(const double *values)
{
    _lambda_optimizer.prepare_calculation(values);
    _gamma_optimizer.prepare_calculation(values + _p_lambda->count());

}

void gamma_lambda_optimizer::report_precalculation()
{
    LOG(INFO) << "Attempting lambda: " << *_p_lambda << ", alpha: " << _gamma_optimizer.get_alpha();
}

/// results consists of the desired number of lambdas and one alpha value
void gamma_lambda_optimizer::finalize(double *results) {
    _lambda_optimizer.finalize(results);
    _gamma_optimizer.finalize(results + _p_lambda->count());
}

