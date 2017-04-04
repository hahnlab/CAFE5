#ifndef clade_h
#define clade_h

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

class newick_parser; // forward declaration (so it can see friend)

class clade {

 friend newick_parser; // allows newick_parser to set parameter values
  
 private:
  string taxon_name;
  long branch_length; // or lambda value
  clade *p_parent; // needs to be pointer; instance creates infinite loop
  vector<clade*> descendants; // same as above
  vector<clade*>::iterator desc_it, desc_end; // iterator (and its end check) for descendants
  /* methods */
  void name_interior_clade();

 public:
  /* methods */
  clade(): p_parent(NULL), branch_length(0) {} // basic constructor

 clade(string name, long length): taxon_name(name), branch_length(length) {} // constructor giving taxon name and branch length

  ~clade(); // destructor
  
  clade *get_parent();

  void add_descendant(clade *p_descendant);

  void print_immediate_descendants();

  void print_clade();

  bool is_leaf();

  void add_leaf_names(vector<string>& v);

  int get_branch_length() const { return branch_length;  }

  template <typename func> 
  void apply_to_descendants(func f)
  {
	  std::for_each(descendants.begin(), descendants.end(), f);
  }

  

  // fill_internal_node_name method HERE
  
};

#endif
