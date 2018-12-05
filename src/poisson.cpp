#include <vector>
#include <cmath>
#include <random>

#include "poisson.h"
#include "utils.h"
#include "clade.h"
#include "fminsearch.h"
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

class poisson_scorer : public optimizer_scorer
{
    vector<int> sizes;
public:
    poisson_scorer(const vector<gene_family>& gene_families)
    {
        for (auto &fam : gene_families)
        {
            for (auto species : fam.get_species())
                if (fam.get_species_size(species) > 0)
                    sizes.push_back(fam.get_species_size(species) - 1);
        }
    }

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double my_random = distribution(randomizer_engine);
        return std::vector<double>{my_random};
    }
    virtual double calculate_score(double * values) override
    {
        lnLPoisson(values, &sizes);
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


};

vector<double> find_poisson_lambda(vector<gene_family> gene_families)
{
    poisson_scorer scorer(gene_families);
    optimizer opt(&scorer);
    
    auto result = opt.optimize();

    cout << "Empirical Prior Estimation Result : (" << result.num_iterations << " iterations)" << endl;
    cout << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score << endl;

    return result.values;
}

