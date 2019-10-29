#include <iostream>
#include <algorithm>
#include <regex>

#include "clade.h"
#include "gene_family.h"

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
        for (auto desc : _descendants) {
            desc->add_leaf_names(vector_names);
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

    const clade *p_descendant;
    auto descendant_finder = [some_taxon_name, &p_descendant](const clade *clade) {
        if (clade->get_taxon_name() == some_taxon_name)
            p_descendant = clade;
    };

  apply_prefix_order(descendant_finder);

  return p_descendant;
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
    for (auto name : descendant_names) {
        _taxon_name += name;
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

std::map<std::string, int> clade::get_lambda_index_map()
{
    std::map<std::string, int> node_name_to_lambda_index;
    
    auto fn = [&node_name_to_lambda_index](const clade *c) { 
        node_name_to_lambda_index[c->get_taxon_name()] = c->get_lambda_index() - 1; // -1 is to adjust the index offset
    };

    apply_prefix_order(fn);
	return node_name_to_lambda_index;
}

void clade::write_newick(ostream& ost, std::function<std::string(const clade *c)> textwriter) const
{
    if (is_leaf()) {
        ost << textwriter(this);
    }
    else {
        ost << '(';

        // some nonsense to supress trailing comma
        for (size_t i = 0; i< _descendants.size() - 1; i++) {
            _descendants[i]->write_newick(ost, textwriter);
            ost << ',';
        }

        _descendants[_descendants.size() - 1]->write_newick(ost, textwriter);
        ost << ')' << textwriter(this);
    }
}

string clade_index_or_name(const clade* node, const cladevector& order)
{
    auto id = distance(order.begin(), find(order.begin(), order.end(), node));
    if (node->is_leaf())
        return node->get_taxon_name() + "<" + to_string(id) + ">";
    else
    {
        return "<" + to_string(id) + ">";
    }
}

std::set<double> clade::get_branch_lengths() const
{
    set<double> result;
    auto branch_length_func = [&result](const clade* c) { 
        if (c->get_branch_length() > 0.0)
            result.insert(c->get_branch_length()); 
    };
    apply_prefix_order(branch_length_func);
    return result;
}

void clade::validate_lambda_tree(const clade* p_lambda_tree) const
{
	auto g = [](set<string>& s, const clade* c) {
		s.insert(c->get_taxon_name());
	};
	set<string> my_taxa;
	apply_prefix_order([g, &my_taxa](const clade* c) { g(my_taxa, c);  });

	set<string> lambda_taxa;
	p_lambda_tree->apply_prefix_order([g, &lambda_taxa](const clade* c) { g(lambda_taxa, c);  });

	if (my_taxa != lambda_taxa)
	{
		throw std::runtime_error("The lambda tree structure does not match that of the tree");
	}
}

void clade::apply_to_descendants(std::function<void(const clade*)> f) const {

	// apply f to direct descendants
	// could replace with apply_prefix_order for functions f that recur through descendants
	//for_each(_descendants.begin(), _descendants.end(), f); // for_each from std
	// for_each apparently passes by value
	for (auto desc : _descendants)
		f(desc);
}

//! apply the functor f to this clade and also to all descendants.
void clade::apply_prefix_order(std::function<void(const clade*)> f) const {
	std::stack<const clade*> stack;
	stack.push(this);
	while (!stack.empty())
	{
		auto c = stack.top();
		stack.pop();

		// Moving from right to left in the tree because that's what CAFE does
		auto it = c->_descendants.rbegin();
		for (; it != c->_descendants.rend(); ++it)
		{
			stack.push(*it);
		}
		f(c);
	}
}

//! apply the functor f to this clade and also to all descendants, by starting
// with the leaf nodes and moving up the tree
void clade::apply_reverse_level_order(std::function<void(const clade*)> f) const {
	std::stack<const clade*> stack;
	std::queue<const clade*> q;

	q.push(this);
	while (!q.empty())
	{
		/* Dequeue node and make it current */
		auto current = q.front();
		q.pop();
		stack.push(current);

		for (auto i : current->_descendants)
		{
			/* Enqueue child */
			q.push(i);
		}
	}

	while (!stack.empty())
	{
		auto current = stack.top();
		stack.pop();
		f(current);
	}
}

clade* parse_newick(std::string newick_string, bool parse_to_lambdas) {

	std::regex tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;");

	auto new_clade = [](clade* p_parent) {
		clade* p_new_clade = new clade();
		if (p_parent != NULL) {
			p_new_clade->_p_parent = p_parent;
		}

		return p_new_clade;
	};

	int lp_count(0), rp_count(0);
	sregex_iterator regex_it(newick_string.begin(), newick_string.end(), tokenizer);
	sregex_iterator regex_it_end;
	clade* p_root_clade = new_clade(NULL);
	p_root_clade->is_lambda_clade = parse_to_lambdas; // if user does not provide lambda for root, we need to make the root specifically a lambda clade if we are parsing to lambdas

	clade* p_current_clade = p_root_clade; // current_clade starts as the root

	// The first element below is empty b/c I initialized it in the class body
	for (; regex_it != regex_it_end; regex_it++) {
		/* Checking all regex */
		// cout << regex_it->str() << endl;

		/* Start new clade */
		if (regex_it->str() == "(") {
			/* Checking '(' regex */
			// cout << "Found (: " << regex_it->str() << endl;

			p_current_clade = new_clade(p_current_clade); // move down the tree (towards the present)
			p_current_clade->get_parent()->add_descendant(p_current_clade); // can't forget to add the now current clade to its parent's descendants vector
			lp_count++;
		}

		else if (regex_it->str() == ",") {
			/* Checking ',' regex */
			// cout << "Found ,: " << regex_it->str() << endl;

			/* The if block below is for when the newick notation omits the external parentheses, which is legal */
			if (p_current_clade == p_root_clade) {
				cout << "Found root!" << endl;
				p_root_clade = new_clade(NULL);
				p_current_clade->_p_parent = p_root_clade; // note that get_parent() cannot be used here because get_parent() copies the pointer and it would be the copy that would be assigned p_root_clade... and then the copy would just be thrown away
				p_current_clade->get_parent()->add_descendant(p_current_clade);
			}

			/* Start new clade at same level as the current clade */
			p_current_clade = new_clade(p_current_clade->get_parent()); // move to the side of the tree
			p_current_clade->get_parent()->add_descendant(p_current_clade); // adding current clade as descendant of its parent
		}

		/* Finished current clade */
		else if (regex_it->str() == ")") {
			/* checking ')' regex */
			// cout << "Found ): " << regex_it->str() << endl;

			p_current_clade = p_current_clade->get_parent(); // move up the tree (into the past)
			rp_count++;
		}

		/* Finished newick string */
		else if (regex_it->str() == ";") {
			/* Checking ';' regex */
			// cout << "Found ;: " << regex_it->str() << endl;
			break;
		}

		/* Reading branch length */
		else if (regex_it->str()[0] == ':') {
			/* Checking ':' regex */
			// cout << "Found :: " << regex_it->str() << endl;

			if (parse_to_lambdas)
			{
				int ind = strtol(regex_it->str().substr(1).c_str(), nullptr, 0);
				p_current_clade->_lambda_index = ind;
				p_current_clade->is_lambda_clade = true;
			}
			else
			{
				p_current_clade->_branch_length = atof(regex_it->str().substr(1).c_str()); // atof() converts string into float
				p_current_clade->is_lambda_clade = false;
			}
		}

		/* Reading taxon name */
		else {
			/* Checking species name string regex */
			// cout << "Found species name: " << regex_it->str() << endl;

			p_current_clade->_taxon_name = regex_it->str();
			clade* p_parent = p_current_clade->get_parent();
			/* If this species has a parent, we need to update the parent's name */
			if (p_parent != NULL) {
				p_parent->_name_interior_clade(); // update parent's name, _name_interior_clade() is a void method
			}
		}
	}

	// since user is not required to set a lambda index for the root, go ahead and assign it to the first lambda
	// so the rest of the code doesn't get confused
	if (p_root_clade->is_lambda_clade)
	{
		if (p_root_clade->get_lambda_index() == 0)
			p_root_clade->_lambda_index = 1;

		auto validator = [](const clade* c) {
			if (c->_lambda_index < 1)
				throw std::runtime_error("Invalid lambda index set for " + c->get_taxon_name());
		};
		p_root_clade->apply_reverse_level_order(validator);
	}
	return p_root_clade;
}
