#include "child_calculator.h"
#include "lambda.h"

using namespace std;

//! Operator () overload.
/*!
Makes child_calculator a functor.
The functor is called once on each child by likelihood_computer through apply_to_descendants().
Note that depending on whether one or multiple lambdas are specified, the computation of the likelihood will be different. It is the abstract class lambda (which has a pure virtual method calculate_child_factor) that decides how to do it.
*/
void child_calculator::operator()(const clade * child) {
  if (_probabilities[child].empty())
  {
    throw std::runtime_error("Child node probabilities not calculated");
  }

  _factors[child].resize(_probabilities_vec_size);
  // cout << "Child factor size is " << _probabilities_vec_size << endl;
  _factors[child] = _lambda->calculate_child_factor(_calc, child, _probabilities[child], s_min_family_size, s_max_family_size, c_min_family_size, c_max_family_size);

  // p(node=c,child|s) = p(node=c|s)p(child|node=c) integrated over all c
  // remember child likelihood[c]'s never sum up to become 1 because they are likelihoods conditioned on c's.
  // incoming nodes to don't sum to 1. outgoing nodes sum to 1
}

//! Method.
/*!
Called by likelihood_computer after all children have been processed. It multiplies all factors together and updates the _probabilities map.
*/
void child_calculator::update_probabilities(const clade *node) {
  vector<double>& node_probs = _probabilities[node];
  node_probs.resize(_probabilities_vec_size);

  for (int i = 0; i < node_probs.size(); i++) {
    node_probs[i] = 1;
    auto it = _factors.begin();

    for (; it != _factors.end(); it++) {
      node_probs[i] *= it->second[i];
    }
  }
}

