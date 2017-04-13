#ifndef clade_h
#define clade_h

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stack>

/* Ask Ben:
1) branch_length is long, but in get_branch_length()'s declaration we use int
2) Why the const in get_branch_length() and get_taxon_name()
3) How is it that print_immediate_descendants can do (*desc_it)->taxon_name, if taxon_name is private?
4) We need to make branch lengths "long floats". How do we do this?
*/

using namespace std;

/* Forward declaration of newick_parser class, so class clade can see friend */
class newick_parser; // actual declaration in utils.h

/* Main clade class that will hold tree structures */
class clade {

 friend newick_parser; // allows newick_parser to set parameter values
  
 private:
   string _taxon_name;

   double _branch_length; // or lambda value

   clade *_p_parent; // needs to be pointer; instance creates infinite loop

   vector<clade*> _descendants; // same as above

   vector<clade*>::iterator _desc_it, _desc_end; // iterator (and its end check) for descendants

   /* methods */
   void _name_interior_clade();

 public:
   /* methods */
   clade(): _p_parent(NULL), _branch_length(0) {} // basic constructor

   clade(string some_taxon_name, double length): _taxon_name(some_taxon_name), _branch_length(length) {} // constructor giving taxon name and branch length

   ~clade(); // destructor
  
   clade *get_parent();

   void add_descendant(clade *p_descendant);

   void add_leaf_names(vector<string>& vector_names);

   bool is_leaf();

   bool is_root();
   
   double get_branch_length() const { return _branch_length; }

   vector<clade*> find_internal_nodes();
   
   clade *find_descendant(string some_taxon_name);

   double find_branch_length(string some_taxon_name);

   string get_taxon_name() const { return _taxon_name; }

   void print_immediate_descendants(); // for testing purposes as of now

   void print_clade(); // for testing purposes as of now

   template <typename func> void apply_to_descendants(func f) {
     // apply f to direct descendants
     // could replace with apply_prefix_order for functions f that recur through descendants
     for_each(_descendants.begin(), _descendants.end(), f); // for_each from std
   }

   template <typename func> void apply_prefix_order(func& f) { // f must be passed by reference to avoid copies being made of f 
     // having a copy made would mean any state variables of f would be lost
     std::stack<clade *> stack;
     stack.push(this);
     while (!stack.empty())
     {
       clade *c = stack.top();
       stack.pop();

       // Moving from right to left in the tree because that's what CAFE does
       std::vector<clade*>::reverse_iterator it = c->_descendants.rbegin();
       for (; it != c->_descendants.rend(); ++it)
       {
         stack.push(*it);
       }
       f(c);
     }
   }
};

/* This class will store a descendant clade if it finds the provided taxon_name */
class descendant_finder {

 private:
  string _some_taxon_name;
  clade *_p_descendant_found;

 public:
  descendant_finder(string some_taxon_name) : _some_taxon_name(some_taxon_name), _p_descendant_found(NULL) { }
  
  void operator()(clade *clade) {
    if (clade->get_taxon_name() == _some_taxon_name) {
      _p_descendant_found = clade;
    }
  }
  
  clade *get_result() { return _p_descendant_found;  }
};
#endif
