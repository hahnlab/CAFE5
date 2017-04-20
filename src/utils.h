#ifndef utils_h
#define utils_h

#include <regex>
#include <string>
#include <iostream>
#include <fstream>

class clade;

class GeneFamily
{
  std::map<std::string, int> species_size;
public:
  int get_species_size(std::string species) const
  {
    return species_size.at(species);
  }
  void set_species_size(std::string species, int sz)
  {
    species_size[species] = sz;
  }
};

class newick_parser {

 private:
  int lp_count, rp_count;
  std::regex tokenizer;

  /* methods */
  clade *new_clade(clade *p_parent);
  
 public:
  std::string newick_string;

  /* methods */
  newick_parser(): tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;"), lp_count(0), rp_count(0) {} // constructor
  clade *parse_newick();
};

class likelihood_computer
{
  // represents probability of the node having various family sizes
  std::map<clade *, std::vector<double> > _probabilities;
  GeneFamily *_family;
  int _max_possible_family_size;
  double _lambda;
public:
  likelihood_computer(int max_possible_family_size, double lambda, GeneFamily *family) : _max_possible_family_size(max_possible_family_size), _lambda(lambda)
  {
    _family = family;
  }
  void operator()(clade *node);

  double *get_likelihoods() const { return NULL; }
};

using pair_type = map<int, int>::value_type;

bool max_value(const pair_type & p1, const pair_type & p2); // def'n in utils.cpp
#endif
