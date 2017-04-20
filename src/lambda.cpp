#include <vector>
#include <map>
#include <cmath>
#include "clade.h"
#include "fminsearch.h"
#include "utils.h"

using namespace std;

struct lambda_args
{
  clade *tree;
  vector<double> prior_rfsize; // prior is a poisson distribution on the root size based on leaves' size
  vector<GeneFamily> gene_families;
};

double lambda_searcher(double* plambda, void* v)
{
  lambda_args *args = static_cast<lambda_args *>(v);

  int gene_family_count = 0;  // holds number of gene families
  clade *tree = args->tree;

  vector<double> ML(gene_family_count);
  vector<double> MAP(gene_family_count);
  int i, j;
  double score = 0;
  for (i = 0; i < args->gene_families.size(); i++)	// i: family index
  {
    GeneFamily& fam = args->gene_families[i];
    map<clade *, int> node_family_sizes;  // TODO: this holds the family count at each node
    int root_family_size = 0;   // TODO
    int max_possible_family_size = 0;  // TODO

    likelihood_computer pruner(max_possible_family_size, 0.001, &fam);
    tree->apply_reverse_level_order(pruner);
    vector<double> likelihood = pruner.get_likelihoods(tree);		// likelihood of the whole tree = multiplication of likelihood of all nodes

    ML[i] = *std::max_element(likelihood.begin(), likelihood.end());	// this part find root size condition with maxlikelihood for each family			

                                                                          // get posterior by adding lnPrior to lnLikelihood
    vector<double> posterior(root_family_size);
    // prior is a poisson distribution on the root size based on leaves' size
    for (j = 0; j < root_family_size; j++)	// j: root family size
    {
      // likelihood and posterior both starts from 1 instead of 0 
      posterior[j] = exp(log(likelihood[j]) + log(args->prior_rfsize[j]));	//prior_rfsize also starts from 1
    }

    MAP[i] = *std::max_element(posterior.begin(), posterior.end());			// this part find root size condition with maxlikelihood for each family			
    if (ML[i] == 0)
    {
      score = log(0);
      break;
    }
    score += log(MAP[i]);			// add log-posterior across all families
  }
  return -score;
}

void set_prior_rfsize_empirical(vector<double>& rfsize)
{

}

void find_best_lambda(clade *tree)
{
  lambda_args args;
  args.tree = tree;
  set_prior_rfsize_empirical(args.prior_rfsize);
  FMinSearch *pfm = fminsearch_new_with_eq(lambda_searcher, 1, &args);
  pfm->tolx = 1e-6;
  pfm->tolf = 1e-6;
  fminsearch_min(pfm, NULL);
  double *re = fminsearch_get_minX(pfm);
}