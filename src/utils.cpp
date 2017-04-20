#include "clade.h"
#include "utils.h"
#include "fminsearch.h"
#include "probability.h"

clade *newick_parser::parse_newick() {

  sregex_iterator regex_it(newick_string.begin(), newick_string.end(), tokenizer);
  sregex_iterator regex_it_end;
  clade *p_root_clade = new_clade(NULL);
  clade *p_current_clade = p_root_clade; // current_clade starts as the root
  
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
      
      p_current_clade->_branch_length = atof(regex_it->str().substr(1).c_str()); // atof() converts string into float
   }
    
    /* Reading taxon name */
    else {
      /* Checking species name string regex */
      // cout << "Found species name: " << regex_it->str() << endl;

      p_current_clade->_taxon_name = regex_it->str();
      clade *p_parent = p_current_clade->get_parent();
      /* If this species has a parent, we need to update the parent's name */
      if (p_parent != NULL) {
	p_parent->_name_interior_clade(); // update parent's name, _name_interior_clade() is a void method
      }
    }
  }

  return p_root_clade;
}

clade *newick_parser::new_clade(clade *p_parent) {

  clade *p_new_clade = new clade();
  if (p_parent != NULL) {
    p_new_clade->_p_parent = p_parent;
  }

  return p_new_clade;
}

vector<vector<double> > get_matrix(int size, int branch_length, double lambda)
{
  vector<vector<double> > result(size);
  for (int s = 0; s < size; s++)
  {
    result[s].resize(size);
    for (int c = 0; c < size; c++)
    {
      result[s][c] = the_probability_of_going_from_parent_fam_size_to_c(lambda, branch_length, s, c);
      //cout << "s = " << s << " c= " << c << ", result=" << result[s][c] << endl;
    }
  }
  return result;
}

vector<double> matrix_multiply(const vector<vector<double> >& matrix, const vector<double>& v)
{
  vector<double> result(matrix.size());

  for (int s = 0; s < matrix.size(); s++)
  {
    result[s] = 0;
    for (int c = 0; c < matrix[s].size(); c++)
    {
      result[s] += matrix[s][c] * v[c];
    }
  }

  return result;
}

ostream& operator<<(ostream& ost, const vector<vector<double> >& matrix)
{
  for (int s = 0; s < matrix.size(); s++)
  {
    for (int c = 0; c < matrix[s].size(); c++)
    {
      ost << matrix[s][c] << "\t";
    }
    ost << "\n";
  }

  return ost;
}

class child_calculator {
  map<clade *, vector<double> > _factors;
  map<clade *, vector<double> >& _probabilities;
  int _max_root_family_size;
  double _lambda;
public:
  child_calculator(int max_root_family_size, double lambda, map<clade *, vector<double> >& probabilities) : _max_root_family_size(max_root_family_size), _lambda(lambda), _probabilities(probabilities)
  {
      
  }

  int num_factors() { return _factors.size();  }
  void operator()(clade * child)
  {
    cout << "() Calc is " << this << endl;

    _factors[child].resize(_max_root_family_size);
    vector<vector<double> > matrix = get_matrix(_factors[child].size(), child->get_branch_length(), _lambda);
    //cout << "Child matrix calculated for " << child->get_taxon_name() << endl << matrix << endl << "Done" << endl;

    _factors[child] = matrix_multiply(matrix, _probabilities[child]);
    for (int i = 0; i < _factors[child].size(); ++i)
      cout << "  Factor for child " << child->get_taxon_name() << i << ": " << _factors[child][i] << endl;
    // p(node=c,child|s) = p(node=c|s)p(child|node=c) integrated over all c
    // remember child likelihood[c]'s never sum up to become 1 because they are likelihoods conditioned on c's.
    // incoming nodes to don't sum to 1. outgoing nodes sum to 1
    cout << "Factor count: " << _factors.size() << endl;
  }

  void update_probabilities(clade *node)
  {
    cout << "Updating probabilities for " << node->get_taxon_name() << " - " << _factors.size() << " factors" << endl;
    _probabilities[node].resize(_max_root_family_size);
    for (int i = 0; i < _probabilities[node].size(); i++)
    {
      _probabilities[node][i] = 1;
      map<clade *, std::vector<double> >::iterator it = _factors.begin();
      for (; it != _factors.end(); it++)
      {
        cout << "    Multiplying by " << it->second[i];
        _probabilities[node][i] *= it->second[i];
      }
    }
  }
};

void likelihood_computer::operator()(clade *node)
{
  cout << "Calculating probabilities for " << node->get_taxon_name() << endl;
  if (node->is_leaf())
  {
    _probabilities[node].resize(_max_possible_family_size);
    int species_size = _family->get_species_size(node->get_taxon_name());

    _probabilities[node][species_size] = 1.0;

  }
  else
  {
    child_calculator calc(_max_possible_family_size, _lambda, _probabilities);
    cout << "l Calc is " << &calc << endl;
    node->apply_to_descendants(calc);
    cout << "l Calc is " << &calc << endl;
    cout << "Factors: " << calc.num_factors() << endl;
    calc.update_probabilities(node);
    cout << "Factors: " << calc.num_factors() << endl;
  }
  cout << "Probabilities calculated for " << node->get_taxon_name() << endl;
}

bool max_value(const pair_type & p1, const pair_type & p2) { 

  return p1.first < p2.first; 
}
