#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>

/// Base class for use by the optimizer
class optimizer_scorer {
public:
    virtual std::vector<double> initial_guesses() = 0;

    virtual double calculate_score(double *values) = 0;

    virtual void finalize(double *result) = 0;

    virtual ~optimizer_scorer() {}
};

class error_model;
class lambda;
class model;
class root_equilibrium_distribution;
class clade;
class base_model;

/// 
class lambda_optimizer : public optimizer_scorer
{
    const clade *_p_tree;
    lambda *_p_lambda;
    base_model *_p_model;
    root_equilibrium_distribution *_p_distribution;
public:
    lambda_optimizer(const clade *p_tree, lambda *p_lambda, base_model* p_model, root_equilibrium_distribution *p_distribution) :
        _p_tree(p_tree),
        _p_lambda(p_lambda),
        _p_model(p_model),
        _p_distribution(p_distribution),
        quiet(false)
    {
    }

    std::vector<double> initial_guesses();

    virtual double calculate_score(double *values);

    virtual void finalize(double *results);

    bool quiet;
};


/// optimize lambdas and epsilons together
class lambda_epsilon_optimizer : public optimizer_scorer
{
    error_model* _p_error_model;
    lambda *_p_lambda;
    model *_p_model;
    root_equilibrium_distribution *_p_distribution;

    double _longest_branch;
    std::vector<double> current_guesses;
public:
    lambda_epsilon_optimizer(model* p_model,
        error_model *p_error_model,
        root_equilibrium_distribution* p_distribution,
        lambda *p_lambda,
        double longest_branch) :
        _p_model(p_model),
        _p_error_model(p_error_model),
        _p_distribution(p_distribution),
        _p_lambda(p_lambda),
        _longest_branch(longest_branch)
    {
#ifdef SILENT
        quiet = true;
#endif
    }

    std::vector<double> initial_guesses();

    virtual double calculate_score(double *values);
    virtual void finalize(double *results);

    bool quiet = false;
};

class gamma_model;

class gamma_lambda_optimizer : public optimizer_scorer
{
    gamma_model *_p_model;
    root_equilibrium_distribution *_p_distribution;
    const clade *_p_tree;
    lambda *_p_lambda;
public:
    gamma_lambda_optimizer(const clade *p_tree, lambda *p_lambda, gamma_model * p_model, root_equilibrium_distribution *p_distribution) :
        _p_tree(p_tree),
        _p_lambda(p_lambda),
        _p_model(p_model),
        _p_distribution(p_distribution)
    {

    }

    std::vector<double> initial_guesses();

    double calculate_score(double *values);

    /// results consists of the desired number of lambdas and one alpha value
    void finalize(double *results);
};

class gamma_optimizer : public optimizer_scorer {
    gamma_model *_p_model;
    root_equilibrium_distribution* _p_distribution;
public:
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual double calculate_score(double * values) override;
    virtual void finalize(double * result) override;
    gamma_optimizer(gamma_model* p_model, root_equilibrium_distribution* prior) : _p_model(p_model), _p_distribution(prior)
    {

    }
};




#endif
