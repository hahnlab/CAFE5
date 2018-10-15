#include "clade.h"
#include "gene_family.h"
#include "utils.h"
#include "io.h"

using namespace std;

/* Recursive destructor */
clade::~clade() {

  for (auto i : _descendants) {
    delete i; // recursively delete all descendants
  }
}

clade *clade::get_parent() const {

  return _p_parent; // memory address
}

double clade::get_branch_length() const
{
    if (is_lambda_clade)
        throw std::runtime_error("Requested branch length from lambda tree");

    return _branch_length;
}

int clade::get_lambda_index() const
{
    if (!is_lambda_clade)
        throw std::runtime_error("Requested lambda index from branch length tree");

    return _lambda_index;
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
vector<const clade*> clade::find_internal_nodes() const {

  vector<const clade*> internal_nodes;

  /* Base case: returns empty vector */
  if (is_leaf()) { return internal_nodes; }

  else {
    internal_nodes.push_back(this);

    for (auto i : _descendants) {
      auto descendant = i->find_internal_nodes(); // recursion
      if (!descendant.empty()) { internal_nodes.insert(internal_nodes.end(), descendant.begin(), descendant.end()); }
    }

    return internal_nodes;
  }
}


  /* Recursively find pointer to clade with provided taxon name */
const clade *clade::find_descendant(string some_taxon_name) const {
  descendant_finder finder(some_taxon_name);
  apply_prefix_order(finder);
  return finder.get_result();
}
  
/* Finds branch length of clade with provided taxon name. Does so by calling find_descendant(), which recursively searches the tree */
double clade::find_branch_length(string some_taxon_name) {

  auto clade = find_descendant(some_taxon_name);
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

class named_lambda_index_getter
{
public:
	std::map<std::string, int> node_name_to_lambda_index;
	void operator()(const clade *c)
	{
		node_name_to_lambda_index[c->get_taxon_name()] = c->get_lambda_index()-1; // -1 is to adjust the index offset
	}
};

std::map<std::string, int> clade::get_lambda_index_map()
{
	named_lambda_index_getter m;
	apply_prefix_order(m);
	return m.node_name_to_lambda_index;
}

void clade::write_newick(ostream& ost, std::function<std::string(const clade *c)> textwriter) const
{
    if (is_leaf()) {
        ost << textwriter(this);
    }
    else {
        ost << '(';

        // some nonsense to supress trailing comma
        for (int i = 0; i< _descendants.size() - 1; i++) {
            _descendants[i]->write_newick(ost, textwriter);
            ost << ',';
        }

        _descendants[_descendants.size() - 1]->write_newick(ost, textwriter);
        ost << ')' << textwriter(this);
    }
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
