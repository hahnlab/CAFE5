#ifndef FAMILY_H
#define FAMILY_H

#include <vector>
#include <iosfwd>
#include <map>
//#include "probability.h"

class clade;
class lambda;
class matrix_cache;
class error_model;

typedef std::map<clade *, int> trial;

trial * simulate_family_from_root_size(clade *tree, int root_family_size, int max_family_size, lambda * p_lambda, error_model *p_error_model);

// std::vector<std::vector<trial *> > simulate_families_from_distribution(clade *p_tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda, error_model *p_error_model);

std::vector<trial *> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, lambda *p_lambda, error_model *p_error_model);

class random_familysize_setter {
private:
    trial *_p_tth_trial;
    int _max_family_size; //!< We simulate from 0 to _max_family_size-1
    lambda *_p_lambda;
    matrix_cache* _calculator; //!< Does the birth-death model computations
    error_model *_p_error_model;

public:
    //! Constructor
    random_familysize_setter(trial *p_tth_trial, int max_family_size, lambda * p_lambda, matrix_cache* p_calc, error_model *error) :
        _p_tth_trial(p_tth_trial), _max_family_size(max_family_size), _p_lambda(p_lambda), _calculator(p_calc), _p_error_model(error) {
    }

    //! Operator () overload.
    /*!
    Makes random_family_setter a functor.
    The functor is called once on each node of the tree, recursively, and in each recursion it populates _p_tth_trial.
    */
    void operator()(clade *node);

};


#endif
