#include <vector>
#include <map>
#include <random>
#include "clade.h"
#include "probability.h"
#include "family_generator.h"
#include "matrix_cache.h"

extern std::mt19937 randomizer_engine;

//! Set the family size of a node to a random value, using parent's family size
/*!
  Starting from 0, the gene family size of the child (c) is increased until the cumulative probability of c (given the gene family size s of the parent) exceeds a random draw from a uniform distribution. When this happens, the last c becomes the child's gene family size.
 
  Note that the smaller the draw from the uniform, the higher the chance that c will be far away from s.

  The birth-death probability calculations are done using a class (probability_calculator defined in probability.h.
*/
void random_familysize_setter::operator()(const clade *node) {

    if (node->is_root()) { return; } // if node is root, we do nothing

                                   /* Drawing random number from uniform */
                                   //  std::default_random_engine gen(static_cast<long unsigned int>(time(0)));
                                   // double rnd = dis(gen);
    int parent_family_size = (*_p_tth_trial)[node->get_parent()];
    size_t c = 0; // c is the family size we will go to

    double lambda = _p_lambda->get_value_for_clade(node);
    double branch_length = node->get_branch_length();

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    c = select_size(parent_family_size, lambda, branch_length, distribution(randomizer_engine));

    if (node->is_leaf() && _p_error_model != NULL)
    {
        if (c >= _p_error_model->get_max_count())
        {
            throw runtime_error("Trying to simulate leaf family size that was not included in error model");
        }
        auto probs = _p_error_model->get_probs(c);
       
        double rnd = distribution(randomizer_engine);
        if (rnd < probs[0])
        {
            c--;
        }
        else if (rnd > (1 - probs[2]))
        {
            c++;
        }
    }

    (*_p_tth_trial)[node] = c;
}

/// returns the index of the item for which the cumulative sum of probabilities of a smaller size is less than rnd
int random_familysize_setter::select_size(int parent_family_size, double lambda, double branch_length, double rnd)
{
    int c = 0;
    double cumul = 0;
    if (parent_family_size > 0) {
        auto m = _cache.get_matrix(branch_length, lambda);
        if (m.is_zero_except_00())  // saturated, return an invalid value
            return parent_family_size;

        for (; c < _max_family_size - 1; c++) { // Ben: why -1
            double prob = m.get(parent_family_size, c);
            //double prob = the_probability_of_going_from_parent_fam_size_to_c(lambda, branch_length, parent_family_size, c);
            cumul += prob;

            if (cumul >= rnd) {
                break;
            }
        }
    }
    return c;
}

