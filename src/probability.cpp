#include <algorithm>
#include <cmath>
#include <iostream>
#include "utils.h" // for gene_family class
#include "clade.h"
#include "probability.h"

/* Useful links
1) http://www.rskey.org/gamma.htm # explanation for lgamma
2) http://www.physics.unlv.edu/~pang/cp_c.html # c code
*/

/* Necessary for old C implementation of what now is lgamma
#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
*/

/* Necessary for old C implementation of what now is lgamma 
double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
 24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
 -5.395239384953e-6 };
*/

/* Old C implementation of what now is lgamma */
/*
double gammaln(double a)
{
	int n;
	double p = __Qs[0];
	double a_add_5p5 = a + 5.5;
	for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
	return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}
*/


/* START: Math tools --------------------- */
/* Old C implementation necessary for set_node_familysize_random. Now using uniform_real_distribution()
*/
double unifrnd()
{
  double result = rand() / (RAND_MAX + 1.0); // rand() returns an int from 0 to RAND_MAX (which is defined in std); the +1.0 is there probably so that we do not draw exactly 1.
  return result;
}

double chooseln(double n, double r)
{

  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);
  return lgamma(n + 1) - lgamma(r + 1) - lgamma(n - r + 1);
}

/* END: Math tools ----------------------- */

/* START: Birth-death model components --------------------- */

/* Eqn. (1) in 2005 paper. Assumes u = lambda */
double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff)
{

  int m = std::min(c, s);
  double lastterm = 1;
  double p = 0.0;
  int s_add_c = s + c;
  int s_add_c_sub_1 = s_add_c - 1;
  int s_sub_1 = s - 1;

  for (int j = 0; j <= m; j++) {
    double t = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
    p += (exp(t) * lastterm); // Note that t is in log scale, therefore we need to do exp(t) to match Eqn. (1)
    lastterm *= coeff; // equivalent of ^j in Eqn. (1)
  }

  return std::max(std::min(p, 1.0), 0.0);
}

double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size)
{
  //std::cout << "Lambda: " << lambda << ", branch_length: " << branch_length;
  double alpha = lambda*branch_length / (1 + lambda*branch_length);
  double coeff = 1 - 2 * alpha;
  
  //printf("Birthdeath rate for 1, 0 (alpha=%f, coeff=%f), : %f\n", alpha, coeff, birthdeath_rate_with_log_alpha(1, 0, log(alpha), coeff));

  double result = birthdeath_rate_with_log_alpha(parent_size, size, log(alpha), coeff);

//  if (result < .000000000000000001)
//  {
//    std::cout << "result= " << result << " lambda=" << lambda << " branch_length:" << branch_length << " From " << parent_size << " to " << size << std::endl;
//  }
  return result;
}

/* END: Birth-death model components ----------------------- */

/* START: Likelihood computation ---------------------- */

//! Compute transition probability matrix for all gene family sizes from 0 to size-1 (=_max_root_family_size-1)
vector<vector<double> > get_matrix(int size, int branch_length, double lambda) {
    
    cout << "Computing matrix for branch " << branch_length << " lambda " << lambda << endl;

    vector<vector<double> > result(size);
    for (int s = 0; s < size; s++) {
        result[s].resize(size);
    
        for (int c = 0; c < size; c++) {
            result[s][c] = the_probability_of_going_from_parent_fam_size_to_c(lambda, branch_length, s, c);
            //cout << "s = " << s << " c= " << c << ", result=" << result[s][c] << endl;
        }
    }
    
    return result;
}

//! Take in a matrix and a vector, compute product, return it
vector<double> matrix_multiply(const vector<vector<double> >& matrix, const vector<double>& v) {
    vector<double> result(matrix.size());

    for (int s = 0; s < matrix.size(); s++) {
        result[s] = 0;
    
        for (int c = 0; c < matrix[s].size(); c++) {
            result[s] += matrix[s][c] * v[c];
        }
    }

    return result;
}

//! Computation and storage of the vector of likelihoods (one likelihood/family size) of an internal node
/*!
  This class computes the vector of probabilities (of the data given the model; i.e., the likelihoods) for an internal node.
   
  An instance of it is created by likelihood_computer for a focal internal node (in the code: "node").
 
  likelihood_computer then calls child_calculator as a function (using the () operator overload) on each descendant of the internal node (via apply_to_descendants()), computing the vector of likelihoods for each one. This vector (called a "factor") is then stored in the _factors map as keys, at _factors[child]. There will be as many factors as there are children.
  
  likelihood_computer then calls the update_probabilities() method of child_calculator, which multiplies the factors together, and stores the result (a vector of likelihoods) in _probabilities[node] ("node" here is the focal internal node).
 
  Note that the _probabilities map stores the vectors of likelihoods from all nodes in the tree, and is a member of the likelihood_computer class.
*/
class child_calculator {
private:
    map<clade *, vector<double> > _factors; //!< keys = pointers to clade objects (children of internal node), values = clade's contribution (factor) to the vector of likelihoods of the internal node
    map<clade *, vector<double> >& _probabilities; //!< (member of likelihood_computer) keys = pointer to clade object (all nodes in the tree), values = clade's vector of likelihoods
    int _max_root_family_size; //!< max gene family size for which likelihood is to be computed
    double _lambda; //!< mambda used in likelihood computation

public:
    //! Constructor.
    /*!
      Used once per internal node by likelihood_computer().
    */
    child_calculator(int max_root_family_size, double lambda, map<clade *, vector<double> >& probabilities) : _max_root_family_size(max_root_family_size), _lambda(lambda), _probabilities(probabilities) {
    }

    int num_factors() { return _factors.size(); }
  
    //! Operator () overload.
    /*!
      Makes child_calculator a functor. 
      The functor is called once on each child by likelihood_computer through apply_to_descendants().
     */
    void operator()(clade * child) {
        _factors[child].resize(_max_root_family_size); // Ben: Maybe I'm wrong, but are we ignoring the largest possible gene family size?
        // Ben: Because if you look at get_matrix(), s goes from 0 to size-1. So if _max_root_family_size is 10, our vector stops at 9, yes?
        vector<vector<double> > matrix = get_matrix(_factors[child].size(), child->get_branch_length(), _lambda); // Ben: is _factors[child].size() the same as _max_root_family_size? If so, why not use _max_root_family_size instead?

        _factors[child] = matrix_multiply(matrix, _probabilities[child]);
    // p(node=c,child|s) = p(node=c|s)p(child|node=c) integrated over all c
    // remember child likelihood[c]'s never sum up to become 1 because they are likelihoods conditioned on c's.
    // incoming nodes to don't sum to 1. outgoing nodes sum to 1
    }

    //! Method.
    /*!
     Called by likelihood_computer after all children have been processed. It multiplies all factors together and updates the _probabilities map.
    */
    void update_probabilities(clade *node) {
        _probabilities[node].resize(_max_root_family_size);
    
        for (int i = 0; i < _probabilities[node].size(); i++) {
            _probabilities[node][i] = 1;
            map<clade *, std::vector<double> >::iterator it = _factors.begin();
            
            for (; it != _factors.end(); it++) {
                _probabilities[node][i] *= it->second[i];
            }
        }
    }
};

//! Operator () overload for likelihood_computer.
/*!
  The operator () overload here allows likelihood_computer to be called as a function (i.e., it makes likelihood_computer a functor).
  This is what allows likelihood_computer to be called recursively in the pruning algorithm, through apply_reverse_level_order.
*/
void likelihood_computer::operator()(clade *node) {
    if (node->is_leaf()) {
        _probabilities[node].resize(_max_possible_family_size);
        int species_size = _family->get_species_size(node->get_taxon_name());
        _probabilities[node][species_size] = 1.0;
    }
  
    else {
        child_calculator calc(_max_possible_family_size, _lambda, _probabilities);
        node->apply_to_descendants(calc);
        calc.update_probabilities(node);
    }
}

/* END: Likelihood computation ---------------------- */

