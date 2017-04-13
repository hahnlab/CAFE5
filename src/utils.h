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
public:
  likelihood_computer(int max_possible_family_size, GeneFamily *family) : _max_possible_family_size(max_possible_family_size)
  {
    _family = family;
  }
  void operator()(clade *node);

  double *get_likelihoods() const { return NULL; }
};


#endif
