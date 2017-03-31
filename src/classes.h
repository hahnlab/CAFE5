#ifndef classes_h
#define classes_h

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "prefix", optional_argument, NULL, 'p' },
  { "suffix", optional_argument, NULL, 's' },
  { "outfmt", optional_argument, NULL, 'o' },
  { 0, 0, 0, 0 }
};

class clade {

 public:
  string taxon_name;
  long branch_length; // or lambda value
  clade *p_parent; // needs to be pointer; instance creates infinite loop
  vector<clade*> descendants; // same as above
  vector<clade*>::iterator desc_it, desc_end; // iterator (and its end check) for descendants

  /* methods */
  
  clade(): p_parent(NULL), branch_length(0) {} // basic constructor
  clade(string name, long length): taxon_name(name), branch_length(length) {} // constructor giving taxon name and branch length
  ~clade(); // destructor
  
  clade *get_parent();

  void add_descendant(clade *p_descendant);

  void print_immediate_descendants();

  void print_clade();

  void name_interior_clade();

  bool am_leaf();

  void add_leaf_names(vector<string>& v);

  // fill_internal_node_name method HERE
  
};

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

void clade::add_leaf_names(vector<string>& v)
{
	if (descendants.empty())
		v.push_back(taxon_name);
	else
	{
		for (int i = 0; i < descendants.size(); ++i)
		{
			descendants[i]->add_leaf_names(v);
		}
	}
}

void clade::name_interior_clade() {
	vector<string> desc_names;
	add_leaf_names(desc_names);
	std::sort(desc_names.begin(), desc_names.end());

	taxon_name.clear();
	for (int i = 0; i < desc_names.size(); ++i)
		taxon_name += desc_names[i];

	cout << "Interior node renamed: " << taxon_name << endl;
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

bool clade::am_leaf() {

  if (descendants.empty()) { return true; }
  else { return false; }
}
#endif
