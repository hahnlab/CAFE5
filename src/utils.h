#ifndef utils_h
#define utils_h

#include <regex>
#include <string>
#include <iostream>
#include <fstream>

string tree_str = "((Sp_A:1.0,Sp_B:1.0):1,Sp_C:2);";
regex tokenizer("\\(|\\)|[^\\s\\(\\)\\[\\]\\\'\\:\\;\\,]+|\\:[+-]?[::digit::]*\\.?[::digit::]+([eE][+-]?[::digit::]+)?|\\,|\\[(\\\\.|[^\\]])*\\]|\\\'(\\\\.|[^\\\'])*\\\'|\\;");

sregex_iterator regex_it(tree_str.begin(), tree_str.end(), tokenizer);
sregex_iterator regex_it_end; // an uninitialized sregex_iterator represents the ending position!

class newick_parser {

 public:
  string newick_string;
  int lp_count, rp_count;

 /* methods */
  clade *parse_newick();
  clade *new_clade(clade *p_parent);
};


clade *newick_parser::parse_newick() {

  clade *p_root_clade = new_clade(NULL);
  clade *p_current_clade = p_root_clade;
  
  // The first element below is empty b/c I initialized it in the class body
  for (; regex_it != regex_it_end; regex_it++) {
    
  /* cout << regex_it->str() << endl; */
    if (regex_it->str() == "(") {
      // cout << regex_it->str() << endl;
      p_current_clade = new_clade(p_current_clade);
      lp_count++;
    }

    else if (regex_it->str() == ",") {
      // cout << regex_it->str() << endl;
      /* if current clade is root */
      if (p_current_clade == p_root_clade) {
	p_root_clade = new_clade(NULL);
	p_current_clade->p_parent = p_root_clade;
      }

      /* start new clade at same level as the current clade */
      /* p_parent = process_clade(p_current_clade); */
    }
    
  }

  return NULL;
}

clade *newick_parser::new_clade(clade *p_parent) {

  clade *p_new_clade = new clade();
  if (p_parent != NULL) {
    p_new_clade->p_parent = p_parent;
  }

  return p_new_clade;
}

/* clade *newick_parser::process_clade(clade *some_clade) { */
  
/* } */
#endif
