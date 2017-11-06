#include <getopt.h>
#include <cmath>
#include <map>

#include "io.h"
#include "execute.h"
#include "utils.h"
#include "clade.h"
#include "probability.h"
#include "family_generator.h"
#include "fminsearch.h"
#include "poisson.h"
#include "core.h"
#include "gamma.h"
#include "fminsearch.h"
#include "lambda.h"
#include "gamma_core.h"
#include "root_equilibrium_distribution.h"
#include "base_model.h"

using namespace std;

vector<double> get_posterior(vector<gene_family> gene_families, int max_family_size, int max_root_family_size, double lambda, clade *p_tree)
{
  vector<double> posterior(max_family_size);

  vector<double> root_poisson_lambda = find_poisson_lambda(gene_families);

  vector<double> prior_rfsize = get_prior_rfsize_poisson_lambda(0, max_family_size, root_poisson_lambda[0]);
  //for (int i = 0; i < prior_rfsize.size(); ++i)
  //  cout << "prior_rfsize " << i << "=" << prior_rfsize[i] << endl;

  probability_calculator calc;
  single_lambda lam(&calc, lambda);
  likelihood_computer pruner(max_root_family_size, max_family_size, &lam, &gene_families[0]);

  p_tree->apply_reverse_level_order(pruner);
  //cout << "Pruner complete" << endl;
  vector<double> likelihood = pruner.get_likelihoods(p_tree);		// likelihood of the whole tree = multiplication of likelihood of all nodes

  int root_max_family_size = 30;
  for (int i = 0; i < root_max_family_size-1; i++) {
    // our lk vector has a lk for family size = 0, but we are ignoring it (the root cannot have size 1) -- so we do i + 1 when calculating the posterior
    // our prior on the root family size starts at family size = 1, and goes to root_max_family size; because there is no 0 here, we do root_max_family_size - 1 in the for loop
      posterior[i] = exp(log(likelihood[i+1]) + log(prior_rfsize[i]));	//prior_rfsize also starts from 1
	  cout << "Likelihood of size " << i << " is " << likelihood[i + 1] << " - prior " << prior_rfsize[i] << endl;
  }

  return posterior;
}

void do_something(clade *c)
{
	cout << "Clade id is: " << c->get_taxon_name() << endl;
}

double pvalue(double v, vector<double>& conddist)
{
	int idx = std::upper_bound(conddist.begin(), conddist.end(), v) - conddist.begin();
	return  idx / (double)conddist.size();
}

void call_viterbi(int max_family_size, int max_root_family_size, int number_of_simulations, lambda *p_lambda, vector<gene_family>& families, clade* p_tree)
{
	double lambda_val = dynamic_cast<single_lambda *>(p_lambda)->get_single_lambda();
	auto cd = get_conditional_distribution_matrix(p_tree, max_root_family_size, max_family_size, number_of_simulations, lambda_val);

	for (int i = 0; i < families.size(); ++i)
	{
		int max = families[i].get_max_size();
		max_root_family_size = rint(max*1.25);

		likelihood_computer pruner(max_root_family_size, max_family_size, p_lambda, &families[i]);
		p_tree->apply_reverse_level_order(pruner);
		auto lh = pruner.get_likelihoods(p_tree);
		std::cout << "likelihoods = ";
		for (int i = 0; i < lh.size(); ++i)
			std::cout << lh[i] << " ";
		std::cout << std::endl;

		double observed_max_likelihood = pruner.max_likelihood(p_tree);	// max value but do we need a posteriori value instead?
	//	std::cout << "observed_likelihood = ";
	//	for (int i = 0; i < root_family_size; ++i)
	//		std::cout << observed_likelihood[i] << " ";
	//	std::cout << std::endl;
	//	max_family_size = std::max(50, max_family_size / 5);
		vector<double> pvalues(max_root_family_size);
		for (int s = 0; s < max_root_family_size; s++)
		{
			pvalues[s] = pvalue(observed_max_likelihood, cd[s]);
		}
	}
}

int cafexp(int argc, char *const argv[]) {
    /* START: Option variables for main() */
    int args; // getopt_long returns int or char
    int prev_arg;
    srand(10);
    int max_family_size = -1; // needed for defining matrix size
    int max_root_family_size = -1; // needed for defining matrix size
    input_parameters my_input_parameters;
    
    clade *p_tree = new clade();
    std::vector<gene_family> gene_families;
    probability_calculator calculator; // for computing lks
    map<int, int>* p_rootdist_map = NULL; // for sims

    while (prev_arg = optind, (args = getopt_long(argc, argv, "i:o:t:y:n:f:l:m:k:a:s::g::p::", longopts, NULL)) != -1 ) {
    // while ((args = getopt_long(argc, argv, "i:t:y:n:f:l:e::s::", longopts, NULL)) != -1) {
        if (optind == prev_arg + 2 && *optarg == '-') {
            cout << "You specified option " << argv[prev_arg] << " but it requires an argument. Exiting..." << endl;
            exit(EXIT_FAILURE);
            // args = ':';
            // --optind;
        }
            
        switch (args) {
            case 'i':
                my_input_parameters.input_file_path = optarg;
                break;
            case 'o':
                my_input_parameters.output_prefix = optarg;
                break;
            case 't':
		my_input_parameters.tree_file_path = optarg;
                break;
            case 'y':
		my_input_parameters.lambda_tree_file_path = optarg;
                break;
            case 's':
                // Number of fams simulated defaults to 0 if -f is not provided
                if (optarg != NULL) { my_input_parameters.nsims = atoi(optarg); }
                break;
            case 'l':
		my_input_parameters.fixed_lambda = atof(optarg);
                break;
            case 'p':
                my_input_parameters.use_uniform_eq_freq = false; // If the user types '-p', the root eq freq dist will not be a uniform
                // If the user provides an argument to -p, then we do not estimate it
                if (optarg != NULL) { my_input_parameters.poisson_lambda = atof(optarg); }
                break;
            case 'm':
		my_input_parameters.fixed_multiple_lambdas = optarg;
                break;
            case 'k':
                if (optarg != NULL) { my_input_parameters.n_gamma_cats = atoi(optarg); }
                break;
            case 'a':
                my_input_parameters.fixed_alpha = atof(optarg);
                break;
//            case 'n':
//		my_input_parameters.nsims = atoi(optarg);
//                break;
            case 'f':
		my_input_parameters.rootdist = optarg;
                break;
            case 'g':
		my_input_parameters.do_log = true;
                break;
            case ':':   // missing argument
                fprintf(stderr, "%s: option `-%c' requires an argument",
                        argv[0], optopt);
                break;
            default: // '?' is parsed (
                cout << "Exiting..." << endl;
                return EXIT_FAILURE; //abort ();
          }
      }
 
    execute my_executer;
    try {
        my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

        /* -t */
        // clade *p_tree = my_executer.read_input_tree(my_input_parameters); // phylogenetic tree
        my_executer.read_input_tree(my_input_parameters, p_tree); // phylogenetic tree

        /* -i */
        if (!my_input_parameters.input_file_path.empty()) {
            // Populates (pointer to) vector of gene family instances, max_family_size and max_root_family_size (last two passed by reference)
            my_executer.read_gene_family_data(my_input_parameters, max_family_size, max_root_family_size, p_tree, &gene_families);
        }
            
        /* -y */
        clade *p_lambda_tree = my_executer.read_lambda_tree(my_input_parameters);

        /* -l/-m (in the absence of -l, estimate) */
        lambda *p_lambda = my_executer.read_lambda(my_input_parameters, calculator, p_lambda_tree);

        root_equilibrium_distribution* p_prior = root_eq_dist_factory(my_input_parameters, &gene_families);
        
        // When computing or simulating, only base or gamma model is used. When estimating, base and gamma model are used (to do: compare base and gamma w/ LRT)
        // Build model takes care of -f
        vector<model *> models = build_models(my_input_parameters, p_tree, p_lambda, &gene_families, max_family_size, max_root_family_size);
        
        // -f cannot have been specified
        if (my_input_parameters.rootdist.empty()) {
        
            // If lambda was fixed, compute!
            if (p_lambda) {            
                my_executer.compute(models, &gene_families, p_prior, my_input_parameters, max_family_size, max_root_family_size); // passing my_input_parameters as reference
            }
        
            // If lambda was not fixed, estimate!
            else {
                for (model* p_model : models) {
                    p_lambda = my_executer.estimate_lambda(p_model, my_input_parameters, p_prior, p_tree, p_lambda_tree, &gene_families, max_family_size, max_root_family_size, calculator);
                    p_model->set_lambda(p_lambda);
                }
            
            // Printing: take estimated values, re-compute them for printing purposes
            my_executer.compute(models, &gene_families, p_prior, my_input_parameters, max_family_size, max_root_family_size); 
            }
        }
        delete p_prior;

        /* -s */
        if (my_input_parameters.nsims != 0 || !my_input_parameters.rootdist.empty()) {
            // -s is provided an argument (-f is not), using -i to obtain root eq freq distr'n
            if (my_input_parameters.nsims != 0) {
                // place holder for estimating poisson lambda if -p, or using uniform as root eq freq distr'n
            }

            // -f is provided (-s does not have an argument), not using -i
            else if (!my_input_parameters.rootdist.empty()) {
                cout << "Using -f, not using -i, nsims = " << my_input_parameters.nsims << endl;
                my_executer.simulate(models, my_input_parameters);
            }
        }

        /* -g */
        if (my_input_parameters.do_log) {

            string prob_matrix_suffix = "_tr_prob_matrices.txt";
            string prob_matrix_file_name = my_input_parameters.output_prefix + prob_matrix_suffix;
            std::ofstream ofst(my_input_parameters.output_prefix + prob_matrix_suffix);
            calculator.print_cache(ofst, max_family_size);
        }
        /* END: Printing log file(s) */

        return 0;
    }
    catch (runtime_error& err) {
        cout << err.what() << endl;
        return EXIT_FAILURE;
    }
} // end main

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

  // cout << "Testing pruning algorithm" << endl;
  // GeneFamily family;
  // likelihood_computer pruner(600, &family);
  // p_tree->apply_prefix_order(pruner);
  //double* likelihood = pruner.get_likelihoods();		// likelihood of the whole tree = multiplication of likelihood of all nodes


/* START: Testing implementation of clade class */

/* Initializing tree */
/*
  parser.newick_string = "((A:1,B:1):2,C:3);";
  clade *p_tree = parser.parse_newick();
*/

/* Testing print_immediate_descendants */
/*
cout << "Testing print_immediate_descendants():" << endl;
  p_tree->print_immediate_descendants();
*/

/* Testing print_clade() method */
/*
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
