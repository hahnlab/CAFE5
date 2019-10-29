#ifndef ERROR_MODEL_H
#define ERROR_MODEL_H

#include <vector>
#include <string>
#include <map>

/* START: Reading in error model file
In this file, \b maxcnt is the largest family size observed in the dataset.Errorclasses(for all following rows) are
defined with \b cntdiffand act as labels for error distributions for each gene family size.Error classes must be
space - delimited positive or negative integers(and 0).The error class with label 0 means that this corresponds to no
change in gene family size due to error.After the first two lines, each possible family size in the dataset(size 0 to
    \b maxcnt) should have an error distribution defined.Any omitted family size follows the distribution for the previous
    row.The error distribution for each count should be space delimited probabilities whose columns correspond to the
    error classes defined in line two.

    <em>Default< / em>: No error model is applied.

    \note
    1. You should not specify any negative error correction for family size of 0 as this cannot occur(i.e., there
        can't be negative gene family sizes);
        2. The rows of the error model file must sum to 1;
3. If any gene counts are missing from the error model file, CAFE will assume the same error distribution from the
previous line.This can also be used as a shortcut if you know that all of the gene counts are specified with the
same error distribution : simply enter the first four lines(\b maxcnt, \b cntdiff, <em>family size < / em >= 0, 1) into
    the error model file and CAFE will use the distribution for family size = 1 as the distribution for all gene family sizes.
    */
class error_model {
private:
    size_t _max_family_size;

    std::vector<int> _deviations; //!< Deviations from the true gene family (e.g., -1 0 1)

    std::vector<std::vector<double> > _error_dists; //!< Each vector element will be a gene family size; the vector of doubles inside (e.g., 0.1 0.8 0.1) will be the probs of deviating from the true value

public:
    error_model();

    //! Set max family size for which deviations apply
    void set_max_family_size(size_t max_cnt);

    //! Set deviations
    void set_deviations(std::vector<std::string> deviations);

    //! Set deviation probability vector for a certain family size
    void set_probabilities(size_t fam_size, std::vector<double>);

    //! Get deviation probability vector for a certain family size
    std::vector<double> get_probs(size_t fam_size) const;

    size_t n_deviations() const {
        return _deviations.size();
    }

    size_t get_max_family_size() const {
        return _error_dists.size();
    }
    std::vector<double> get_epsilons() const;
    void replace_epsilons(std::map<double, double>* new_epsilons);
    void update_single_epsilon(double new_epsilon);

    friend void write_error_model_file(std::ostream& ost, error_model& errormodel);
};

#endif
