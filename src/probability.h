#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size);
double chooseln(double n, double k);
double unifrnd();

class likelihood_computer
{
  // represents probability of the node having various family sizes
  std::map<clade *, std::vector<double> > _probabilities;
  gene_family *_family;
  int _max_possible_family_size;
  double _lambda;
public:
  likelihood_computer(int max_possible_family_size, double lambda, gene_family *family) : _max_possible_family_size(max_possible_family_size), _lambda(lambda)
  {
    _family = family;
  }
  void operator()(clade *node);

  std::vector<double> get_likelihoods(clade *node) const
  { 
    return _probabilities.at(node); 
  }
};

#endif
