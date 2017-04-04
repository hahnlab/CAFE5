#include "utils.h"
#include "clade.h"
#include <getopt.h>
#include <cmath>
#include <map>
#include "probability.h"

using namespace std;

void simulate_families(clade *tree, int num_trials, std::vector<int> root_dist, int max_family_size, double lambda);
extern map<clade *, int> family_sizes;

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
	if (k > n) {
		return 0;
	}
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
  parser.newick_string = "((A:1.0,B:1.0):1.0,C:2.0);";
  clade *p_tree = parser.parse_newick();
  
  /* START: Testing implementation of Clade class */

  /* creating (((A,B)AB,C)ABC) */
  // clade *p_child11 = new clade("A", 1); // (leaf) child11 = A
  // clade *p_child12 = new clade("B", 1); // (leaf) child12 = B
  // clade *p_child1 = new clade("AB", 1); // (internal) child1 = AB
  // p_child1->add_descendant(p_child11);
  // p_child1->add_descendant(p_child12);
  // clade *p_child2 = new clade("C", 2); // (leaf) child2 = C
  // clade parent;
  // parent.taxon_name = "ABC"; // (root) parent = ABC
  // parent.add_descendant(p_child1);
  // parent.add_descendant(p_child2);
  
  /* testing print_immediate_descendants */
  cout << "Testing print_immediate_descendants():" << endl;
  // parent.print_immediate_descendants();
  // p_child1->print_immediate_descendants();
  // p_tree->print_immediate_descendants();
  p_tree->print_immediate_descendants();

  /* testing print_clade() method */
  cout << "Testing print_clade():" << endl;
  // parent.print_clade();
  p_tree->print_clade();

  /* testing am_leaf() method */
  cout << "Testing am_leaf():" << endl;
  // if (!parent.am_leaf()) { cout << "I am not leaf\n"; }
  // if (p_child11->am_leaf()) { cout << "I am leaf\n"; }
  if (!p_tree->is_leaf()) { cout << "I am not leaf\n"; }
  if (!p_tree->is_leaf()) { cout << "I am leaf\n"; }
  
  /* END: Testing implementation of Clade class - */
  
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

  std::vector<int> root_dist;
  root_dist.push_back(5);
  simulate_families(p_tree, 10, root_dist, 10, 0.01);

  map <clade *, int>::iterator it = family_sizes.begin();
  for (; it != family_sizes.end(); ++it)
  {
	  cout << it->first->get_taxon_name() << " family size: " << it->second << endl;
  }
  return 0;
}

