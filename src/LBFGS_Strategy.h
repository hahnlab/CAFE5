#ifndef LBFGS_STRATEGY_H
#define LBFGS_STRATEGY_H

#include "../config.h"

#ifdef HAVE_EIGEN_CORE

#include "optimizer.h"
#include <string>
#include <vector>

#include <Eigen/Core>
#include "LBFGS.h"

class gradient_calculator
{
public:
    candidate _last_success;
    gradient_calculator() : _last_success(0)
    {

    }
    void approximate_infinite_gradient(const std::vector<double>& t, Eigen::VectorXd& grad) const;
    void update(const std::vector<double>& v, double score)
    {
        if (!std::isinf(score))
        {
            _last_success.values = v;
            _last_success.score = score;
        }
    }
    void compute_gradients(const std::vector<double>& t, Eigen::VectorXd& grad, double score, std::function <double(const std::vector<double>& values)> score_func);
};

class LBFGS_strategy : public OptimizerStrategy {
    std::vector<double> _scale_factors;
    gradient_calculator _gradient_calculator;
public:
    virtual void Run(FMinSearch* pfm, optimizer::result& r, std::vector<double>& initial) override;
    virtual std::string Description() const override { return "Broyden-Fletcher-Goldfarb-Shanno algorithm"; };
    double calculate_score_and_gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad, optimizer_scorer* scorer, LBFGSpp::LBFGSSolver<double>& solver);
    double calculate_score_from_scaled_values(optimizer_scorer* scorer, const std::vector<double>& scaled_values);
    void calculate_scale_values(const std::vector<double>& v);

};
#endif

#endif
