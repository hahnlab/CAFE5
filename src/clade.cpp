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
  if (p_parent)
	  p_parent->name_interior_clade();
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

void clade::name_interior_clade() {

  vector<string> desc_names; // vector of names
  add_leaf_names(desc_names); // fills vector of names
  std::sort(desc_names.begin(), desc_names.end()); // sorts alphabetically
  taxon_name.clear(); // resets whatever taxon_name was
  for (int i = 0; i < desc_names.size(); ++i) {
    taxon_name += desc_names[i];
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
