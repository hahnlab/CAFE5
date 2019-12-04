#include <vector>
#include <cmath>
#include <random>

#include "poisson.h"
#include "clade.h"
#include "optimizer.h"
#include "probability.h"
#include "gene_family.h"
#include "optimizer_scorer.h"

extern std::mt19937 randomizer_engine;

using namespace std;

/***********************************************************************
* Poisson Distribution
***********************************************************************/
double poisspdf(int x, double lambda)
{
  return exp(x*log(lambda) - lgamma(x + 1) - lambda);
}

vector<double> get_prior_rfsize_poisson_lambda(int min_family_size, int max_family_size, double poisson_lambda)
{
  int num_sizes = max_family_size - min_family_size;
  vector<double> prior_rfsize(num_sizes);

  for (int i = 0; i<num_sizes; i++) {

    //param->prior_rfsize[i] = poisspdf(param->pcafe->rootfamilysizes[0]+i, parameters[0]);					// poisson
    prior_rfsize[i] = poisspdf(i, poisson_lambda);					// shifted poisson
                                                   //param->prior_rfsize[i] = gampdf(param->pcafe->rootfamilysizes[0]+i, parameters[0], parameters[1]);	// gamma
  }
  return prior_rfsize;
}

poisson_scorer::poisson_scorer(const vector<gene_family>& gene_families)
{
    for (auto &fam : gene_families)
    {
        for (auto species : fam.get_species())
            if (fam.get_species_size(species) > 0)
                leaf_family_sizes.push_back(fam.get_species_size(species) - 1);
    }
}

// Inherited via optimizer_scorer
std::vector<double> poisson_scorer::initial_guesses()
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double my_random = distribution(randomizer_engine);
    return std::vector<double>{my_random};
}

double poisson_scorer::calculate_score(const double * values)
{
    return lnLPoisson(values);
}

double poisson_scorer::lnLPoisson(const double* plambda)
{
    double lambda = plambda[0];
    double score = accumulate(leaf_family_sizes.begin(), leaf_family_sizes.end(), 0.0, [lambda](double x, int sz) {
        double ll = poisspdf((double)sz, lambda);
        if (std::isnan(ll) || std::isinf(ll) || ll == 0) {
            return x;
        }
        return x + log(ll);
        });

    return -score;
}
