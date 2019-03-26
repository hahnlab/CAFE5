#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include "optimizer_scorer.h"

class clade;
class gene_family;

class poisson_scorer : public optimizer_scorer
{
    std::vector<int> leaf_family_sizes;
public:
    poisson_scorer(const std::vector<gene_family>& gene_families);

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual double calculate_score(const double * values) override;

    double lnLPoisson(const double* plambda);
};

std::vector<double> get_prior_rfsize_poisson_lambda(int min_family_size, int max_family_size, double poisson_lambda);

#endif
