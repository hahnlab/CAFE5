#include <vector>
#include <cmath>

#include "poisson.h"
#include "utils.h"
#include "clade.h"
#include "fminsearch.h"
#include "probability.h"

using namespace std;

/***********************************************************************
* Poisson Distribution
***********************************************************************/
double poisspdf(int x, double lambda)
{
  return exp(x*log(lambda) - lgamma(x + 1) - lambda);
}

double lnLPoisson(double* plambda, void* data)
{
  int i = 0;
  double score = 0;
  double lambda = plambda[0];
  std::vector<int> *p_leaf_family_sizes = (std::vector<int> *)data;
  for (i = 0; i<p_leaf_family_sizes->size(); i++) {
    int x = p_leaf_family_sizes->at(i);
    double ll = poisspdf((double)x, lambda);
    if (std::isnan(ll)) {
      ll = 0;
    }
    score += log(ll);
  }

  return -score;
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

vector<double> find_poisson_lambda(vector<gene_family> gene_families)
{
    vector<int> sizes;
    for (auto &fam : gene_families)
    {
        for (auto species : fam.get_species())
            if (fam.get_species_size(species) > 0)
                sizes.push_back(fam.get_species_size(species) - 1);
    }

  FMinSearch* pfm;
  vector<double> result(1);
  pfm = fminsearch_new_with_eq(lnLPoisson, result.size(), &sizes);
  pfm->tolx = 1e-6;
  pfm->tolf = 1e-6;

  double my_random = unifrnd();

  fminsearch_min(pfm, &my_random);

  result[0] = fminsearch_get_minX(pfm)[0];
  cout << "Empirical Prior Estimation Result : (" << pfm->iters << " iterations)" << endl;
  cout << "Poisson lambda: " << result[0] << " &  Score: " << *pfm->fv << endl;

  return result;
}

