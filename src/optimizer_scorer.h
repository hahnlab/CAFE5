#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>
#include <map>

class error_model;
class lambda;
class model;
class root_equilibrium_distribution;
class clade;
class base_model;

/// @brief Base class for use by the optimizer
//! \ingroup optimizer
class optimizer_scorer {
public:
    virtual ~optimizer_scorer() {}

    virtual std::vector<double> initial_guesses() = 0;

    virtual double calculate_score(const double *values) = 0;
};

//! @brief  Scorer that holds a model and calls its inference method
//! for scoring
//! \ingroup optimizer
class inference_optimizer_scorer : public optimizer_scorer {
protected:
    virtual void prepare_calculation(const double *values) = 0;
    virtual void report_precalculation() = 0;

    lambda *_p_lambda;
    model *_p_model;
    root_equilibrium_distribution *_p_distribution;

public:
    inference_optimizer_scorer(lambda *p_lambda, model* p_model, root_equilibrium_distribution *p_distribution) :
        _p_lambda(p_lambda),
        _p_model(p_model),
        _p_distribution(p_distribution),
        quiet(false)
    {
#ifdef SILENT
        quiet = true;
#endif
    }

    virtual ~inference_optimizer_scorer() {}

    double calculate_score(const double *values) ;

    virtual void finalize(double *result) = 0;

    bool quiet;
};

//! @brief Scorer that optimizes for lambda
//! \ingroup optimizer
class lambda_optimizer : public inference_optimizer_scorer
{
    double _longest_branch;

public:
    lambda_optimizer(lambda *p_lambda, model* p_model, root_equilibrium_distribution *p_distribution, double longest_branch) :
        inference_optimizer_scorer(p_lambda, p_model, p_distribution),
        _longest_branch(longest_branch)
    {
    }

    std::vector<double> initial_guesses() override;

    virtual void finalize(double *results) override;

    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;
};


/// @brief Scorer that optimizes lambdas and epsilons together
//! \ingroup optimizer
class lambda_epsilon_optimizer : public inference_optimizer_scorer
{
    lambda_optimizer _lambda_optimizer;
    error_model* _p_error_model;
    std::vector<double> current_guesses;
public:
    lambda_epsilon_optimizer(
        model* p_model,
        error_model *p_error_model,
        root_equilibrium_distribution* p_distribution,
        const std::map<int, int>& root_distribution_map,
        lambda *p_lambda,
        double longest_branch) :
        inference_optimizer_scorer(p_lambda, p_model, p_distribution),
        _lambda_optimizer(p_lambda, p_model, p_distribution, longest_branch),
        _p_error_model(p_error_model)
    {
    }

    std::vector<double> initial_guesses() override;

    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;

    virtual void finalize(double *results) override;
};

class gamma_model;

//! @brief Scorer that optimizes for alpha
//! \ingroup optimizer
//! \ingroup gamma
class gamma_optimizer : public inference_optimizer_scorer {
    gamma_model *_p_gamma_model;
public:
    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual void finalize(double * result) override;
    gamma_optimizer(gamma_model* p_model, root_equilibrium_distribution* prior);

    double get_alpha() const;
};

//! @brief Scorer that optimizes for both lambda and alpha
//! \ingroup optimizer
//! \ingroup gamma
class gamma_lambda_optimizer : public inference_optimizer_scorer
{
    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;
    lambda_optimizer _lambda_optimizer;
    gamma_optimizer _gamma_optimizer;
public:
    gamma_lambda_optimizer(lambda *p_lambda, gamma_model * p_model, root_equilibrium_distribution *p_distribution, double longest_branch);

    std::vector<double> initial_guesses() override;

    /// results consists of the desired number of lambdas and one alpha value
    void finalize(double *results) override;
};




#endif
