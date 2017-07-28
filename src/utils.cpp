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
      
		if (parse_to_lambdas)
			p_current_clade->_lambda_index = atoi(regex_it->str().substr(1).c_str()); // atoi() converts string into int
		else
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

//! Used when simulating gene families (-s)
std::vector<int> vectorize_map(map<int, int> *p_root_dist) {
    std::vector<int> vectorized_map;
    
    for (map<int, int>::iterator it = p_root_dist->begin(); it != p_root_dist->end(); ++it) {
        for (int i=0; i < it->second; ++i) {
            vectorized_map.push_back(it->first);
        }
    }
    
    return vectorized_map;
}
