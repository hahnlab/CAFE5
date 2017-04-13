#include "utils.h"
#include "clade.h"
#include <getopt.h>
#include <cmath>
#include <map>
#include "probability.h"
#include "family_generator.h"

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "prefix", optional_argument, NULL, 'p' },
  { "suffix", optional_argument, NULL, 's' },
  { "outfmt", optional_argument, NULL, 'o' },
  { 0, 0, 0, 0 }
};

string input_file_path;
string prefix;
string suffix;
string outfmt;

unsigned long long choose(unsigned long long n, unsigned long long k) {

  if (k > n) { return 0; }

  unsigned long long r = 1;

  for (unsigned long long d = 1; d <= k; ++d) {
    r *= n--;
    r /= d;
  }

  return r;
}

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
  parser.newick_string = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93)";

  /* START: Testing implementation of clade class */  
  parser.newick_string = "((A:1,B:1):2,C:3);";
  clade *p_tree = parser.parse_newick();
  
  /* END: Testing implementation of clade class - */
  
  //simulate_families();

  int num_trials = 1;
  int root_family_size = 300;
  int max_family_size = 600;
  double lambda = 0.0017;

  cout << "Cafe says this value should be 0.083" << endl;
  cout << the_probability_of_going_from_parent_fam_size_to_c(0.5, .42, 40, 42) << endl;

  cout << "Cafe says this value should be something" << endl;
  cout << the_probability_of_going_from_parent_fam_size_to_c(lambda, .42, 300, 295) << endl;

  //trial simulation = simulate_families_from_root_size(p_tree, num_trials, root_family_size, max_family_size, lambda);

  //print_simulation(simulation, cout);

//   map<int, int> root_size;
// #if 1
//   root_size[300] = 1;
// #else
//   root_size[1] = 6640;
//   root_size[2] = 1641;
//   root_size[3] = 738;
//   root_size[4] = 333;
//   root_size[5] = 172;
//   root_size[6] = 136;
//   root_size[7] = 84;
//   root_size[8] = 54;
//   root_size[9] = 29;
//   root_size[10] = 26;
//   root_size[11] = 19;
//   root_size[12] = 15;
//   root_size[13] = 17;
//   root_size[14] = 9;
//   root_size[15] = 6;
//   root_size[16] = 7;
//   root_size[17] = 5;
//   root_size[18] = 6;
//   root_size[19] = 7;
//   root_size[20] = 4;
//   root_size[21] = 5;
//   root_size[22] = 5;
//   root_size[23] = 3;
//   root_size[24] = 6;
//   root_size[25] = 1;
//   root_size[26] = 4;
//   root_size[27] = 2;
//   root_size[29] = 1;
//   root_size[31] = 2;
//   root_size[32] = 1;
//   root_size[33] = 1;
//   root_size[34] = 1;
//   root_size[35] = 1;
//   root_size[36] = 1;
//   root_size[37] = 1;
//   root_size[39] = 1;
//   root_size[41] = 1;
//   root_size[48] = 1;
//   root_size[49] = 1;
//   root_size[118] = 1;
//   root_size[309] = 1;
// #endif
//   vector<trial> simulations = simulate_families_from_distribution(p_tree, num_trials, root_size, max_family_size, lambda);

//   print_simulation(simulations, cout);

  cout << "Testing pruning algorithm" << endl;
  GeneFamily family;
  likelihood_computer pruner(600, &family);
  p_tree->apply_prefix_order(pruner);
  //double* likelihood = pruner.get_likelihoods();		// likelihood of the whole tree = multiplication of likelihood of all nodes

  return 0;
}

/* START: Testing implementation of clade class */

/* Testing print_immediate_descendants */
/*
cout << "Testing print_immediate_descendants():" << endl;
  p_tree->print_immediate_descendants();

  // Testing print_clade() method
  cout << "Testing print_clade():" << endl;
  p_tree->print_clade();
*/

/* Testing find_internal_node() method */
/*
  vector<clade *> internal_nodes = p_tree->find_internal_nodes();
  cout << "Internal nodes are: ";
  for (vector<clade *>::iterator int_node_it = internal_nodes.begin(); int_node_it != internal_nodes.end(); ++int_node_it)
	cout << (*int_nodeI-t)->get_taxon_name() << ", ";
  cout << endl;
*/

/* Testing find_branch_length() method */
/*
  cout << "Testing find_branch_length():" << endl;

  long branch_length_A = p_tree->find_branch_length("A");
  cout << "The length of A is " << branch_length_A << endl;

  long branch_length_AB = p_tree->find_branch_length("AB");
  cout << "The length of AB is " << branch_length_AB << endl;

  long branch_length_C = p_tree->find_branch_length("C");
  cout << "The length of C is " << branch_length_C << endl;

  long branch_length_ABC = p_tree->find_branch_length("ABC");
  cout << "The length of root ABC is " << branch_length_ABC << endl;

  long branch_length_Z = p_tree->find_branch_length("Z");
  cout << "The length of inexistent Z is " << branch_length_Z << endl;
/*

/* Testing im_leaf() method */
/*
  cout << "Testing am_leaf():" << endl;

  if (!p_tree->is_leaf()) { cout << p_tree->get_taxon_name() << ": I am not leaf\n"; }
  if (p_tree->find_descendant("A")->is_leaf()) { cout << p_tree->find_descendant("A")->get_taxon_name() << ": I am leaf\n"; }
*/

/* END: Testing implementation of clade class - */
