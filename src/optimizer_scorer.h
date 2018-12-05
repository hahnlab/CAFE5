#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>

class error_model;
class lambda;
class model;
class root_equilibrium_distribution;
class clade;
class base_model;

/// Base class for use by the optimizer
class optimizer_scorer {
public:
    virtual ~optimizer_scorer() {}

    virtual std::vector<double> initial_guesses() = 0;

    virtual double calculate_score(double *values) = 0;
};

class inference_optimizer_scorer : public optimizer_scorer {
protected:
    virtual void prepare_calculation(double *values) = 0;
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

    double calculate_score(double *values) ;

    virtual void finalize(double *result) = 0;

    bool quiet;
};

/// 
class lambda_optimizer : public inference_optimizer_scorer
{
    const clade *_p_tree;
    virtual void prepare_calculation(double *values) override;
    virtual void report_precalculation() override;

public:
    lambda_optimizer(const clade *p_tree, lambda *p_lambda, model* p_model, root_equilibrium_distribution *p_distribution) :
        inference_optimizer_scorer(p_lambda, p_model, p_distribution),
        _p_tree(p_tree)
    {
    }

    std::vector<double> initial_guesses();

    virtual void finalize(double *results);

};


/// optimize lambdas and epsilons together
class lambda_epsilon_optimizer : public inference_optimizer_scorer
{
    error_model* _p_error_model;
    double _longest_branch;
    std::vector<double> current_guesses;
public:
    lambda_epsilon_optimizer(model* p_model,
        error_model *p_error_model,
        root_equilibrium_distribution* p_distribution,
        lambda *p_lambda,
        double longest_branch) :
        inference_optimizer_scorer(p_lambda, p_model, p_distribution),
        _p_error_model(p_error_model),
        _longest_branch(longest_branch)
    {
    }

    std::vector<double> initial_guesses();

    virtual void prepare_calculation(double *values) override;
    virtual void report_precalculation() override;

    virtual void finalize(double *results);
};

class gamma_model;

class gamma_lambda_optimizer : public inference_optimizer_scorer
{
    gamma_model *_p_gamma_model;
    const clade *_p_tree;
    virtual void prepare_calculation(double *values) override;
    virtual void report_precalculation() override;
public:
    gamma_lambda_optimizer(const clade *p_tree, lambda *p_lambda, gamma_model * p_model, root_equilibrium_distribution *p_distribution);

    std::vector<double> initial_guesses();

    /// results consists of the desired number of lambdas and one alpha value
    void finalize(double *results);
};

class gamma_optimizer : public inference_optimizer_scorer {
    gamma_model *_p_gamma_model;
public:
    virtual void prepare_calculation(double *values) override;
    virtual void report_precalculation() override;

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual void finalize(double * result) override;
    gamma_optimizer(gamma_model* p_model, root_equilibrium_distribution* prior);
};




#endif
