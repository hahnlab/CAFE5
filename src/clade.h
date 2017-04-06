#ifndef clade_h
#define clade_h

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

/* Ask Ben:
1) branch_length is long, but in get_branch_length()'s declaration we use int
2) Why the const in get_branch_length() and get_taxon_name()
3) How is it that print_immediate_descendants can do (*desc_it)->taxon_name, if taxon_name is private?
4) We need to make branch lengths "long floats". How do we do this?
*/

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

   clade(string some_taxon_name, long length): taxon_name(some_taxon_name), branch_length(length) {} // constructor giving taxon name and branch length

   ~clade(); // destructor
  
   clade *get_parent();

   void add_descendant(clade *p_descendant);

   void add_leaf_names(vector<string>& vector_names);

   bool is_leaf();

   int get_branch_length() const { return branch_length; }

   long find_branch_length(string some_taxon_name);

   string get_taxon_name() const { return taxon_name; }

   void print_immediate_descendants(); // for testing purposes as of now

   void print_clade(); // for testing purposes as of now

   template <typename func> void apply_to_descendants(func f) {
     for_each(descendants.begin(), descendants.end(), f); // for_each from std
   }
};

#endif
