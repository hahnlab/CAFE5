#ifndef utils_h
#define utils_h

#include <regex>
#include <string>
#include <iostream>
#include <fstream>

/* string tree_str = "((Sp_A:1.0,Sp_B:1.0):1,Sp_C:2);"; */
/* regex tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;"); // this cannot go inside the newick parser... ask Ben. */

class newick_parser {

 public:
  string newick_string;
  int lp_count, rp_count;
  regex tokenizer;

 /* methods */
  newick_parser(): tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;"), lp_count(0), rp_count(0) {} // constructor
  clade *parse_newick();
  clade *new_clade(clade *p_parent);
};

clade *newick_parser::parse_newick() {

  sregex_iterator regex_it(newick_string.begin(), newick_string.end(), tokenizer);
  sregex_iterator regex_it_end;
  clade *p_root_clade = new_clade(NULL);
  clade *p_current_clade = p_root_clade; // current_clade starts as the root
  
  // The first element below is empty b/c I initialized it in the class body
  for (; regex_it != regex_it_end; regex_it++) {
    /* cout << regex_it->str() << endl; */

    /* Start new clade */
    if (regex_it->str() == "(") {
      cout << "Found (: " << regex_it->str() << endl;
      p_current_clade = new_clade(p_current_clade); // move down the tree (towards the present)
      // note that calling new_clade(some_clade) returns a new clade with some_clade as its parent
      p_current_clade->get_parent()->add_descendant(p_current_clade); // linking clades must be done both ways (parent -> child here; child <- parent above)
      p_current_clade->get_parent()->print_immediate_descendants();
      lp_count++;
    }

    else if (regex_it->str() == ",") {
      cout << "Found ,: " << regex_it->str() << endl;
      /* if current clade is root */
      if (p_current_clade == p_root_clade) {
      	cout << "Found root!" << endl;
      	p_root_clade = new_clade(NULL);
      	p_current_clade->p_parent = p_root_clade;
      }

      /* start new clade at same level as the current clade */
      p_current_clade = new_clade(p_current_clade->get_parent()); // move to the side of the tree
      p_current_clade->get_parent()->add_descendant(p_current_clade); // adding current clade as descendant of its parent
      p_current_clade->get_parent()->print_immediate_descendants();     
    }

    /* Finished current clade */
    else if (regex_it->str() == ")") {
      cout << "Found ): " << regex_it->str() << endl;
      p_current_clade = p_current_clade->get_parent(); // move up the tree (into the past)
      rp_count++;
    }
    
    /* Finished newick string */
    else if (regex_it->str() == ";") {
      cout << "Found ;: " << regex_it->str() << endl;
      break;
    }

    /* Reading branch length */
    else if (regex_it->str()[0] == ':') {
      cout << "Found :: " << regex_it->str() << endl;
      p_current_clade->branch_length = atof(regex_it->str().substr(1).c_str()); // atof() converts string into float
    }
    
    /* Reading taxon name */
    else {
      cout << "Found species name: " << regex_it->str() << endl;
      p_current_clade->taxon_name = regex_it->str();
	  clade *p_parent = p_current_clade->get_parent();
	  if (p_parent)
	  { 
		  p_parent->name_interior_clade();
		  p_parent->print_immediate_descendants();
	  }
    }
  }

  return p_root_clade;
}

clade *newick_parser::new_clade(clade *p_parent) {

  clade *p_new_clade = new clade();
  if (p_parent != NULL) {
    p_new_clade->p_parent = p_parent;
  }

  return p_new_clade;
}

#endif
