#include "clade.h"

/* Recursive destructor */
clade::~clade() {

  for (desc_it = descendants.begin(), desc_end = descendants.end(); desc_it != desc_end; desc_it++) {
    delete *desc_it; // desc_it is a pointer, *desc_it is a clade; this statement calls the destructor and makes it recursive
  }
}

/* Returns pointer to parent */
clade *clade::get_parent() {

  return p_parent; // memory address
}

/* Adds descendant to vector of descendants */
void clade::add_descendant(clade *p_descendant) {
	
  descendants.push_back(p_descendant);
  name_interior_clade();
  if (p_parent != NULL) {
    p_parent->name_interior_clade();
  }
}

/* Recursively fills vector of names provided as argument */
void clade::add_leaf_names(vector<string> &vector_names) {

  if (descendants.empty()) {
    vector_names.push_back(taxon_name); // base case (leaf), and it starts returning
  }

  else {
    for (int i = 0; i < descendants.size(); ++i) {
	descendants[i]->add_leaf_names(vector_names);
    }
  }
}

/* Recursively finds branch length of some_taxon_name */
long clade::find_branch_length(string some_taxon_name) {

  long some_branch_length;

  /* Base case: found some_taxon name and is not root */
  if ((taxon_name == some_taxon_name) && (p_parent != NULL)) {
    return branch_length;
  }
  
  /* If reached wrong leaf */
  else if (descendants.empty()) { return 0; }

  else {
    for (desc_it = descendants.begin(), desc_end = descendants.end(); desc_it != desc_end ; desc_it++) {
      some_branch_length = (*desc_it)->find_branch_length(some_taxon_name);
      if (some_branch_length == 0) { continue; }
      else { return some_branch_length; }
    }
    
    return 0;
  }
}

void clade::name_interior_clade() {

  vector<string> descendant_names; // vector of names
  add_leaf_names(descendant_names); // fills vector of names
  sort(descendant_names.begin(), descendant_names.end()); // sorts alphabetically (from std)
  taxon_name.clear(); // resets whatever taxon_name was
  for (int i = 0; i < descendant_names.size(); ++i) {
    taxon_name += descendant_names[i];
  }
}

/* Prints names of immediate descendants */
void clade::print_immediate_descendants() {

  cout << "Me: " << taxon_name << " | Descendants: ";
  
  for (desc_it = descendants.begin(), desc_end = descendants.end(); desc_it != desc_end; desc_it++) {
    cout << (*desc_it)->taxon_name << " ";
  }
  
  cout << endl;
}

/* Recursively prints clade */
void clade::print_clade() {

  cout << "My name is: " << taxon_name << endl;

  /* Base case: it is a leaf */
  if (descendants.empty()) {
    return;
  }
  
  for (desc_it = descendants.begin(), desc_end = descendants.end(); desc_it != desc_end ; desc_it++) {
    (*desc_it)->print_clade();
  }
}

bool clade::is_leaf() {

  if (descendants.empty()) { return true; }
  else { return false; }
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
