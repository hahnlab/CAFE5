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
  parser.newick_string = "((A:1,B:1):2,C:3);";
  clade *p_tree = parser.parse_newick();
  
  /* START: Testing implementation of clade class */  
  /* Testing print_immediate_descendants */
  cout << "Testing print_immediate_descendants():" << endl;
  p_tree->print_immediate_descendants();

  /* Testing print_clade() method */
  cout << "Testing print_clade():" << endl;
  p_tree->print_clade();

  /* Testing im_leaf() method */
  cout << "Testing am_leaf():" << endl;
  if (!p_tree->is_leaf()) { cout << "I am not leaf\n"; }
  if (!p_tree->is_leaf()) { cout << "I am leaf\n"; }

  /* Testing find_branch_length() method */
  cout << "Testing find_branch_length():" << endl;
  long branch_length_A = p_tree->find_branch_length("A");
  long branch_length_AB = p_tree->find_branch_length("AB");
  long branch_length_C = p_tree->find_branch_length("C");
  long branch_length_ABC = p_tree->find_branch_length("ABC"); // this is INCORRECT (ask Ben: need to guard against the root)
  long branch_length_Z = p_tree->find_branch_length("Z");
  cout << "The length of A is " << branch_length_A << endl;
  cout << "The length of AB is " << branch_length_AB << endl;
  cout << "The length of C is " << branch_length_C << endl;
  cout << "The length of root ABC is " << branch_length_AB << endl;
  cout << "The length of inexistent Z is " << branch_length_Z << endl;
  
  /* END: Testing implementation of clade class - */
  
  //simulate_families();

  cout << "lgamma(.5) = " << lgamma(.5) << endl;
  cout << "lgamma(1.5) = " << lgamma(1.5) << endl;
  cout << "lgamma(10.5) = " << lgamma(10.5) << endl;
  cout << "lgamma(100.5) = " << lgamma(100.5) << endl;
  cout << "lgamma(1000.5) = " << lgamma(1000.5) << endl;

  cout << "chooseln(10,5) = " << chooseln(10, 5) << endl;
  cout << "chooseln(10,8) = " << chooseln(10, 8) << endl;
  cout << "chooseln(15,5) = " << chooseln(15, 5) << endl;

  cout << "ln(choose(10, 5)) = " << log(choose(10, 5)) << endl;
//  cout << "choose(10, 5) = " << log(choose(10, 5)) << endl;

  cout << "the probability of going from parent 10 to child 5 = " << the_probability_of_going_from_parent_fam_size_to_c(0.01, 1, 10, 5) << endl;
  cout << "the probability of going from parent 10 to child 8 = " << the_probability_of_going_from_parent_fam_size_to_c(0.01, 1, 10, 8) << endl;
  cout << "the probability of going from parent 10 to child 10 = " << the_probability_of_going_from_parent_fam_size_to_c(0.01, 1, 10, 10) << endl;

  int num_trials = 10;
  int root_family_size = 5;
  int max_family_size = 10;
  double lambda = 0.01;

  trial simulation = simulate_families_from_root_size(p_tree, num_trials, root_family_size, max_family_size, lambda);

  print_simulation(simulation, cout);

  vector<int> root_size;
  root_size.push_back(5);
  root_size.push_back(15);
  root_size.push_back(25);
  vector<trial> simulations = simulate_families_from_distribution(p_tree, num_trials, root_size, max_family_size, lambda);

  print_simulation(simulations, cout);

  return 0;
}
