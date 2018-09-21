#ifndef FAMILY_H
#define FAMILY_H

#include <vector>
#include <iosfwd>
#include <map>

class clade;
class lambda;
class matrix_cache;
class error_model;

typedef std::map<const clade *, int> trial;

trial * simulate_family_from_root_size(const clade *tree, int root_family_size, int max_family_size, const lambda * p_lambda, error_model *p_error_model);

std::vector<trial *> simulate_families_from_root_size(const clade *tree, int num_trials, int root_family_size, int max_family_size, const lambda *p_lambda, error_model *p_error_model);

class random_familysize_setter {
private:
    trial *_p_tth_trial;
    int _max_family_size; //!< We simulate from 0 to _max_family_size-1
    const lambda *_p_lambda;
    const error_model *_p_error_model;

public:
    //! Constructor
    random_familysize_setter(trial *p_tth_trial, int max_family_size, const lambda * p_lambda, const error_model *error) :
        _p_tth_trial(p_tth_trial), _max_family_size(max_family_size), _p_lambda(p_lambda), _p_error_model(error) {
    }

    //! Operator () overload.
    /*!
    Makes random_family_setter a functor.
    The functor is called once on each node of the tree, recursively, and in each recursion it populates _p_tth_trial.
    */
    void operator()(const clade *node);

};


#endif
