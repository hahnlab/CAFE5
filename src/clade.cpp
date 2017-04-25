#include "clade.h"
#include <queue>
#include "utils.h"

/* Recursive destructor */
clade::~clade() {

  for (_desc_it = _descendants.begin(), _desc_end = _descendants.end(); _desc_it != _desc_end; _desc_it++) {
    delete *_desc_it; // desc_it is a pointer, *desc_it is a clade; this statement calls the destructor and makes it recursive
  }
}

clade *clade::get_parent() const {

  return _p_parent; // memory address
}

/* Adds descendant to vector of descendants */
void clade::add_descendant(clade *p_descendant) {
	
  _descendants.push_back(p_descendant);
  _name_interior_clade();
  if (!is_root()) {
    _p_parent->_name_interior_clade();
  }
}

/* Recursively fills vector of names provided as argument */
void clade::add_leaf_names(vector<string> &vector_names) {

  if (_descendants.empty()) {
    vector_names.push_back(_taxon_name); // base case (leaf), and it starts returning
  }

  else {
    for (int i = 0; i < _descendants.size(); ++i) {
	_descendants[i]->add_leaf_names(vector_names);
    }
  }
}

/* Recursively finds internal nodes, and returns vector of clades */
vector<clade*> clade::find_internal_nodes() {

  vector<clade*> internal_nodes;

  /* Base case: returns empty vector */
  if (is_leaf()) { return internal_nodes; }

  else {
    internal_nodes.push_back(this);

    for (_desc_it = _descendants.begin(), _desc_end = _descendants.end(); _desc_it != _desc_end; _desc_it++) {
      vector<clade*> descendant = (*_desc_it)->find_internal_nodes(); // recursion
      if (!descendant.empty()) { internal_nodes.insert(internal_nodes.end(), descendant.begin(), descendant.end()); }
    }

    return internal_nodes;
  }
}


  /* Recursively find pointer to clade with provided taxon name */
clade *clade::find_descendant(string some_taxon_name) {
  descendant_finder finder(some_taxon_name);
  apply_prefix_order(finder);
  return finder.get_result();
#if 0
  cout << "Searching for descendant " << some_taxon_name << " in " << get_taxon_name() << endl;

  /* Base case: found some_taxon name and is not root */
  if (_taxon_name == some_taxon_name) { return this; }

  /* If reached (wrong) leaf */
  else if (is_leaf()) { return NULL; }

  else {
    for (_desc_it = _descendants.begin(), _desc_end = _descendants.end(); _desc_it != _desc_end; _desc_it++) {
    clade *p_descendant_to_find = (*_desc_it)->find_descendant(some_taxon_name); // recursion

    if (p_descendant_to_find != NULL) { return p_descendant_to_find; } // recursion is only manifested if finds provided taxon name

    return NULL; // otherwise returns NULL
    }
  }
#endif
}
  
/* Finds branch length of clade with provided taxon name. Does so by calling find_descendant(), which recursively searches the tree */
double clade::find_branch_length(string some_taxon_name) {

  clade *clade = find_descendant(some_taxon_name);
  if (clade == NULL || clade->is_root()) { return 0; } // guarding against root query

  cout << "Found matching clade" << endl;
  return clade->_branch_length;
}

/* Names interior clades, starting from clade of first method call and going up the tree until root */
void clade::_name_interior_clade() {

  vector<string> descendant_names; // vector of names
  add_leaf_names(descendant_names); // fills vector of names
  sort(descendant_names.begin(), descendant_names.end()); // sorts alphabetically (from std)
  _taxon_name.clear(); // resets whatever taxon_name was
  for (int i = 0; i < descendant_names.size(); ++i) {
    _taxon_name += descendant_names[i];
  }
  
  if (_p_parent)
    _p_parent->_name_interior_clade();
}

/* Prints names of immediate descendants */
void clade::print_immediate_descendants() {

  cout << "Me: " << _taxon_name << " | Descendants: ";
  for (_desc_it = _descendants.begin(), _desc_end = _descendants.end(); _desc_it != _desc_end; _desc_it++) {
    cout << (*_desc_it)->_taxon_name << " ";
  }
  
  cout << endl;
}

/* Recursively prints clade */
void clade::print_clade() {

  int depth = 0;
  clade *p_ancestor = get_parent();
  while (p_ancestor) {
    depth++;
    p_ancestor = p_ancestor->get_parent();
  }
  
  string blanks(depth, ' '); // initializing string with the fill constructor (repeat ' ' depth many times, depth will be some integer), this blanks string will help us indent the printing

  cout << blanks << "My name is: " << _taxon_name << endl;

  /* Base case: it is a leaf */
  if (_descendants.empty()) { return; }
  
  for (_desc_it = _descendants.begin(), _desc_end = _descendants.end(); _desc_it != _desc_end ; _desc_it++) {
    (*_desc_it)->print_clade();
  }
}

bool clade::is_leaf() const {

  return _descendants.empty();
}

bool clade::is_root() const {

  return get_parent() == NULL;
}

/* Function print_clade_name() is used in conjunction with apply_reverse_level_order and apply_prefix order for debugging purposes.
   e.g.,
   cout << "Tree " << p_tree->get_taxon_name() << " in reverse order: " << endl;
   p_tree->apply_reverse_level_order(print_name)   
*/
void print_clade_name(clade *clade) {
  cout << clade->get_taxon_name() << " (length of subtending branch: " << clade->get_branch_length() << ")" << endl;
}

void clade::init_gene_family_sizes(const vector<gene_family>& families)
{
  if (is_leaf())
  {
    vector<string> species_names = families[0].get_species();
    for (int i = 0; i<families.size(); ++i)
    {
      const gene_family& f = families[i];
      this->_gene_family_sizes[f.id()] = f.get_species_size(this->get_taxon_name());
    }
  }
  for (int i = 0; i < _descendants.size(); ++i)
    _descendants[i]->init_gene_family_sizes(families);
}

/* Testing implementation of clade class */

/* creating (((A,B)AB,C)ABC) 
clade *p_child11 = new clade("A", 1); // (leaf) child11 = A
clade *p_child12 = new clade("B", 1); // (leaf) child12 = B
clade *p_child1 = new clade("AB", 1); // (internal) child1 = AB
p_child1->add_descendant(p_child11);
p_child1->add_descendant(p_child12);
clade *p_child2 = new clade("C", 2); // (leaf) child2 = C
clade parent;
parent.taxon_name = "ABC"; // (root) parent = ABC
parent.add_descendant(p_child1);
parent.add_descendant(p_child2);
*/

/* testing print_immediate_descendants()
parent.print_immediate_descendants();
p_child1->print_immediate_descendants();
*/

/* testing print_clade()
parent.print_clade();
*/

/* testing is_leaf()
if (!parent.am_leaf()) { cout << "I am not leaf\n"; }
if (p_child11->am_leaf()) { cout << "I am leaf\n"; }
*/
