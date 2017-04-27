#include <vector>
#include <cmath>

#include "poisson.h"
#include "utils.h"
#include "clade.h"
#include "fminsearch.h"
#include "probability.h"

using namespace std;


class leaf_size_gatherer
{
  vector<gene_family> _gene_families;
public:
  leaf_size_gatherer(vector<gene_family> gene_families) : _gene_families(gene_families)
  {

  }
  void operator()(clade *node)
  {
    if (node->is_leaf())
      for (int i = 0; i<_gene_families.size(); ++i)
        sizes.push_back(node->get_gene_family_size(_gene_families[i].id()));
  }
  vector<int> sizes;
};

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
  //printf("lambda: %f (Poisson) & Score: %f\n", lambda, score);	
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

vector<double> find_poisson_lambda(clade* tree, vector<gene_family> gene_families)
{
  // get all leaf sizes
  leaf_size_gatherer loader(gene_families);
  tree->apply_prefix_order(loader);
  for (vector<int>::iterator it = loader.sizes.begin(); it != loader.sizes.end(); ++it)
    cout << "Leaf size is " << (*it) << endl;

  FMinSearch* pfm;
  vector<double> result(1);
  pfm = fminsearch_new_with_eq(lnLPoisson, result.size(), &loader.sizes);
  pfm->tolx = 1e-6;
  pfm->tolf = 1e-6;

  double my_random = unifrnd();

  fminsearch_min(pfm, &my_random);

  result[0] = fminsearch_get_minX(pfm)[0];
  cout << "Score: " << *pfm->fv << endl;

  return result;
}

