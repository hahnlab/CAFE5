#ifndef CHILD_CALCULATOR_H
#define CHILD_CALCULATOR_H

#include <map>
#include <vector>

class clade;
class lambda;
class matrix_cache;

//! Computation and storage of the vector of likelihoods (one likelihood/family size) of an internal node
/*!
This class computes the vector of probabilities (of the data given the model; i.e., the likelihoods) for an internal node.

An instance of it is created by likelihood_computer for a focal internal node (in the code: "node").

likelihood_computer then calls child_calculator as a function (using the () operator overload) on each descendant of the internal node (via apply_to_descendants()), computing the vector of likelihoods for each one. This vector (called a "factor") is then stored in the _factors map as keys, at _factors[child]. There will be as many factors as there are children. The computation of the factors is different depending on whether one or more lambdas are specified -- this is taken care of by the abstract class lambda (and its pure virtual method calculate_child_factor).

likelihood_computer then calls the update_probabilities() method of child_calculator, which multiplies the factors together, and stores the result (a vector of likelihoods) in _probabilities[node] ("node" here is the focal internal node).

Note that the _probabilities map stores the vectors of likelihoods from all nodes in the tree, and is a member of the likelihood_computer class.
*/
class child_calculator {
private:
    std::map<const clade *, std::vector<double> > _factors; //!< keys = pointers to clade objects (children of internal node), values = clade's contribution (factor) to the vector of likelihoods of the internal node
    std::map<const clade *, std::vector<double> >& _probabilities; //!< (member of likelihood_computer) keys = pointer to clade object (all nodes in the tree), values = clade's vector of likelihoods
    int _probabilities_vec_size; //!< size of vector that will store probabilities (likelihoods)
    const lambda* _lambda; //!< lambda used in likelihood computation
    int s_min_family_size; //!< parent min size (this is an index)
    int s_max_family_size; //!< parent max size (this is an index)
    int c_min_family_size; //!< child min size (this is an index)
    int c_max_family_size; //!< child max size (this is an index)
    const matrix_cache& _calc;

public:
    //! Constructor.
    /*!
    Used once per internal node by likelihood_computer().
    */
    child_calculator(int probabilities_vec_size, const lambda* lambda, const matrix_cache& calc, std::map<const clade *, std::vector<double> >& probabilities, int s_min, int s_max, int c_min, int c_max) : _probabilities_vec_size(probabilities_vec_size), _lambda(lambda), _probabilities(probabilities),
    s_min_family_size(s_min), s_max_family_size(s_max), c_min_family_size(c_min), c_max_family_size(c_max),
    _calc(calc) {}

    int num_factors() { return _factors.size(); }

    //! Operator () overload.
    /*!
    Makes child_calculator a functor.
    The functor is called once on each child by likelihood_computer through apply_to_descendants().
    Note that depending on whether one or multiple lambdas are specified, the computation of the likelihood will be different. It is the abstract class lambda (which has a pure virtual method calculate_child_factor) that decides how to do it.
    */
    void operator()(const clade * child);

    //! Method.
    /*!
    Called by likelihood_computer after all children have been processed. It multiplies all factors together and updates the _probabilities map.
    */
    void update_probabilities(const clade *node);
};



#endif