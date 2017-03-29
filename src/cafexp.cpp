#include "classes.h"
#include "utils.h"
#include <getopt.h>
// #include <iostream>
// #include <string>

using namespace std;

// struct option longopts[] = {
//   { "infile", required_argument, NULL, 'i' },
//   { "prefix", optional_argument, NULL, 'p' },
//   { "suffix", optional_argument, NULL, 's' },
//   { "outfmt", optional_argument, NULL, 'o' },
//   { 0, 0, 0, 0 }
// };

string input_file_path;
string prefix;
string suffix;
string outfmt;

int main(int argc, char *const argv[]) {
  
  int args; // getopt_long returns int or char

  while ((args = getopt_long(argc, argv, "i:opq::", longopts, NULL)) != -1) {
    switch (args) {
    case 'i':
      input_file_path = optarg;
      break;
    case 'p':
      prefix = optarg;
      break;
    case 's':
      suffix = optarg;
      break;
    case 'o':
      prefix = optarg;
      break;
    case ':':   // missing argument
      fprintf(stderr, "%s: option `-%c' requires an argument\n",
	      argv[0], optopt);
      break;
    default:
      abort ();
    }
  }
  
  // ifstream input_file(input_file_path.c_str()); /* the constructor for ifstream takes const char*, not string
  // 						   so we need to use c_str() */
  // string line;
  // while (getline(input_file, line)) {
  //   cout << line << endl;
  // }

  newick_parser parser;
  parser.newick_string = "((Sp_A:1.0,Sp_B:1.0):1.0,Sp_C:2.0);";
  clade *p_tree = parser.parse_newick();
  
  /* START: Testing implementation of Clade class */

  /* creating (((A,B)AB,C)ABC) */
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
  
  /* testing print_immediate_descendants */
  cout << "Testing print_immediate_descendants():" << endl;
  parent.print_immediate_descendants();
  p_child1->print_immediate_descendants();
  p_tree->print_immediate_descendants();

  /* testing print_clade() method */
  cout << "Testing print_clade():" << endl;
  parent.print_clade();

  /* testing am_leaf() method */
  cout << "Testing am_leaf():" << endl;
  if (!parent.am_leaf()) { cout << "I am not leaf\n"; }
  if (p_child11->am_leaf()) { cout << "I am leaf\n"; }
  
  /* END: Testing implementation of Clade class - */
  
  return 0;
}
