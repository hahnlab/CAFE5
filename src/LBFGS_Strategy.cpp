#include "LBFGS_Strategy.h"
#include "optimizer_scorer.h"
#include <algorithm>

#include "../config.h"

#ifdef HAVE_EIGEN_CORE
#include <Eigen/Core>
using namespace LBFGSpp;
using namespace std;

const double INFINITY_APPROXIMATION = 50000000;

//! If we are at a discontinuous point, approximate the gradient by calculating the slope between a known good value 
//! and a very large value
void gradient_calculator::approximate_infinite_gradient(const std::vector<double>& t, Eigen::VectorXd& grad) const
{
    double score = _last_success.score;
    auto slope_func = [score](double current, double last) {
        return (score - INFINITY_APPROXIMATION) / (last - current);
    };

    transform(t.begin(), t.end(), _last_success.values.begin(), grad.data(), slope_func);
}

//! descale the values the optimizer is using back to the actual values we need, 
//! and use the scorer to calculate a value
double LBFGS_strategy::calculate_score_from_scaled_values(optimizer_scorer* scorer, const vector<double>& scaled_values)
{
    vector<double> unscaled(scaled_values.size());
    auto descale_func = [](double scaled_value, double factor) { return scaled_value / factor; };

    transform(scaled_values.begin(), scaled_values.end(), _scale_factors.begin(), unscaled.begin(), descale_func);

    double score = scorer->calculate_score(unscaled.data());
    _gradient_calculator.update(scaled_values, score);

    return score;
}

//! Calculate a gradient by adjusting each dimension by a small amount and finding a score at the point
void gradient_calculator::compute_gradients(const std::vector<double>& t, Eigen::VectorXd& grad, double score, std::function <double(const std::vector<double>& values)> score_func)
{
    if (std::isinf(score))
    {
        approximate_infinite_gradient(t, grad);
    }

    int ndim = grad.size();
    vector<double> gradiations(ndim);
    vector<double> gradiated_scores(ndim);

    for (int i = 0; i < ndim; i++) {

        double adjustment = 1e-2 * fabs(t[i]);

        auto adjusted_t = t;
        adjusted_t[i] += adjustment;

        gradiations[i] = adjusted_t[i];
        gradiated_scores[i] = score_func(adjusted_t);
        while (std::isinf(gradiated_scores[i]))
        {
            // we found a discontinuity looking for a gradient point. Go back closer to our valid point until we find a better choice
            adjustment *= 0.8;
            adjusted_t[i] -= adjustment;
            gradiations[i] = adjusted_t[i];
            gradiated_scores[i] = score_func(adjusted_t);
        }
    }
    for (int i = 0; i < ndim; i++)
    {
        grad[i] = (gradiated_scores[i] - score) / gradiations[i];
    }

}

//! Compute the score at point x, and fill out the gradients
double LBFGS_strategy::calculate_score_and_gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad, optimizer_scorer* scorer, LBFGSSolver<double>& solver)
{
    if (std::isnan(x[0]))
        throw std::invalid_argument("Invalid lambda value");

    // Copy the Eigen vector to a std::vector
    vector<double> t(x.data(), x.data()+x.size());

    double score = calculate_score_from_scaled_values(scorer, t);

    _gradient_calculator.compute_gradients(t, grad, score, [this, scorer](const std::vector<double>& v) { return calculate_score_from_scaled_values(scorer, v); });

    return std::isinf(score) ? INFINITY_APPROXIMATION : score;
}

/// invert the initial values and store the inversions
void LBFGS_strategy::calculate_scale_values(const std::vector<double>& v)
{
    _scale_factors.resize(v.size());
    transform(v.begin(), v.end(), _scale_factors.begin(), [](double val) { return 1 / val; });
}

void LBFGS_strategy::Run(FMinSearch* pfm, optimizer::result& r, std::vector<double>& initial)
{
    LBFGSParam<double> param;
    param.delta = 1e-2;
    param.max_iterations = 25;
    param.epsilon = 1e-3;
    param.m = 7;
    param.past = 0;
    LBFGSSolver<double> solver(param);
    optimizer_scorer* scorer = pfm->scorer;
    auto wrapper = [this, scorer, &solver](const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
        return calculate_score_and_gradient(x, grad, scorer, solver);
    };

    calculate_scale_values(initial);

    /// Pass a vector of 1's as our initial scaled values. Our calculated scale_factors will be used 
    /// when computing an actual score
    Eigen::VectorXd x = Eigen::VectorXd::Constant(initial.size(), 1);
    r.num_iterations = solver.minimize(wrapper, x, r.score);

    /// r.score and r.num_iterations are now set.
    /// we just need to descale the values and convert back from the Eigen::vector
    r.values.resize(initial.size());
    transform(_scale_factors.begin(), _scale_factors.end(), x.data(), r.values.begin(), [](double scale_factor, double computed_value) {
        return computed_value / scale_factor;
        });

    // release memory in a way that will make CPPUnit happy
    vector<double>().swap(_gradient_calculator._last_success.values);
    vector<double>().swap(_scale_factors);
}


#endif

