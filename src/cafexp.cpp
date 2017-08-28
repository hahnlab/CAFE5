#include <getopt.h>
#include <cmath>
#include <map>
#include <iomanip>

#include "io.h"
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

/* Ask Ben */
/*
1) How does using pair_type bla works?
2) If I remove the declar'n of max_value from utils.h (but leave the def'n in utils.cpp), cafexp does not compile
*/

using namespace std;

//vector<gene_family> initialize_sample_families()
//{
//  vector<gene_family> gene_families;
//
//  gene_families.push_back(gene_family("ENS01"));
//  gene_families.push_back(gene_family("ENS02"));
//  gene_families.push_back(gene_family("ENS03"));
//  gene_families.push_back(gene_family("ENS04"));
//
//  for (int i = 0; i < gene_families.size(); ++i)
//  {
//    gene_families[i].set_species_size("A", 5);
//    gene_families[i].set_species_size("B", 10);
//    gene_families[i].set_species_size("C", 2);
//    gene_families[i].set_species_size("D", 6);
//  }
//
//  return gene_families;
//
//}

vector<double> get_posterior(vector<gene_family> gene_families, int max_family_size, int max_root_family_size, double lambda, clade *p_tree)
{
  vector<double> posterior(max_family_size);

  srand(10);

  vector<double> root_poisson_lambda = find_poisson_lambda(p_tree, gene_families);

  vector<double> prior_rfsize = get_prior_rfsize_poisson_lambda(0, max_family_size, root_poisson_lambda[0]);
  //for (int i = 0; i < prior_rfsize.size(); ++i)
  //  cout << "prior_rfsize " << i << "=" << prior_rfsize[i] << endl;

  single_lambda lam(lambda);
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

int main(int argc, char *const argv[]) {
    /* START: Option variables for main() */
    int args; // getopt_long returns int or char
    int prev_arg;
    
    std::string input_file_path;
    vector<gene_family> * p_gene_families = new vector<gene_family>; // storing gene family data
    int max_family_size = -1; // largest gene family size among all gene families
	int max_root_family_size = 30;

    std::string tree_file_path;
    
    std::string lambda_tree_file_path;
    clade *p_lambda_tree = new clade(); // lambda tree
    
    std::string rootdist;

    bool estimate = false;
    
    bool simulate = false;
    int nsims = 0; 
    
    double fixed_lambda = 0.0;
    std::string fixed_multiple_lambdas;
    /* END: Option variables for main() */
  
    /* START: Input variables for inference and simulation */
    // std::vector<gene_family> gene_families;
    /* END: Input variables for inference and simulation */
  
    /* START: Option variables for simulations */
    map<int, int>* p_rootdist_map = NULL;
    //double lambda = 0.0017;
    /* END: Option variables for simulations */

    while (prev_arg = optind, (args = getopt_long(argc, argv, "i:t:y:n:f:l:m:e::s::", longopts, NULL)) != -1 ) {
    // while ((args = getopt_long(argc, argv, "i:t:y:n:f:l:e::s::", longopts, NULL)) != -1) {
        if (optind == prev_arg + 2 && *optarg == '-') {
            cout << "You specified option " << argv[prev_arg] << " but it requires an argument. Exiting..." << endl;
            exit(EXIT_FAILURE);
            // args = ':';
            // --optind;
        }
            
        switch (args) {
            case 'i':
                input_file_path = optarg;
                break;
            case 'e':
                estimate = true;
                break;
            case 't':
                tree_file_path = optarg;
                break;
            case 'y':
                lambda_tree_file_path = optarg;
                break;
            case 's':
                simulate = true;
                break;
            case 'l':
                fixed_lambda = atof(optarg);
                break;
            case 'm':
                fixed_multiple_lambdas = optarg;
                break;
            case 'n':
                nsims = atoi(optarg);
                break;
            case 'f':
                rootdist = optarg;
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


//
//    std::vector<double> lambdas = { 0.0, 0.46881494730996, 0.68825840825707 };
//    std::map<clade *, int> lambda_index_map;
//    std::map<std::string, int> node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
//  
//    lambda_index_map[p_tree] = 0;
//
//    multiple_lambda lambda2(node_name_to_lambda_index, lambdas);
//
 
    try {        
        /* START: Checking conflicting options */
        //! The user cannot specify both -e and -l
        if (estimate && fixed_lambda > 0.0) {
            throw runtime_error("You cannot both estimate (-e) and fix the lambda(s) value(s) (-l). Exiting...");
        }
        
        //! The user cannot specify both -l and -y
        if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
            throw runtime_error("You cannot fix one lambda value (-l) and many lambda values (-m). Exiting...");
        }
        /* END: Checking conflicting options */
        
        /* START: Reading tree (-t) */
        clade *p_tree = read_tree(tree_file_path, false); // phylogenetic tree
        /* END: Reading tree */
        
        /* START: Reading gene family data (-i) */
        if (!input_file_path.empty()) {
            p_gene_families = read_gene_families(input_file_path);
            
            // Iterating over gene families to get max gene family size
            for (std::vector<gene_family>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
                int this_family_max_size = it->get_parsed_max_size();
                if (max_family_size < this_family_max_size)
                    max_family_size = this_family_max_size;
            }

            cout << max_family_size << endl;
        }
        /* END: Reading gene family data */
        
		std::map<std::string, int> node_name_to_lambda_index;
        /* START: Reading lambda tree (-y) */
        if (!lambda_tree_file_path.empty()) {
            p_lambda_tree = read_tree(lambda_tree_file_path, true);

			std::map<clade *, int> lambda_index_map;
			node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
			
			p_lambda_tree->print_clade();
        }
        /* END: Reading lambda tree */
        
		lambda *p_lambda = NULL;
		/* START: Read in user-specified multiple lambdas (-m) */
		if (!fixed_multiple_lambdas.empty()) {
			if (lambda_tree_file_path.empty())
				throw runtime_error("You must specify a lambda tree (-y) if you fix multiple lambda values (-m). Exiting...");

			vector<string> lambdastrings = tokenize_str(fixed_multiple_lambdas, ',');
			vector<double> lambdas(lambdastrings.size());

			transform(lambdastrings.begin(), lambdastrings.end(), lambdas.begin(),
				[](string const& val) {return stod(val);});
			p_lambda = new multiple_lambda(node_name_to_lambda_index, lambdas);
		}
		
		/* START: Computing likelihood of user-specified lambda (-l) */
        if (fixed_lambda > 0.0) {
            cout << "Specified lambda (-l): " << fixed_lambda << ". Computing likelihood..." << endl;
            p_lambda = new single_lambda(fixed_lambda);
            // vector<double> posterior = get_posterior((*p_gene_families), max_family_size, fixed_lambda, p_tree);
            // double map = log(*max_element(posterior.begin(), posterior.end()));
            // cout << "Posterior values found - max log posterior is " << map << endl;
			call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
        }

		if (p_lambda != NULL)
		{
			likelihood_computer pruner(max_root_family_size, max_family_size, p_lambda, &(*p_gene_families)[0]); // likelihood_computer has a pointer to a gene family as a member, that's why &(*p_gene_families)[0]
			p_tree->apply_reverse_level_order(pruner);
			vector<double> likelihood = pruner.get_likelihoods(p_tree);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			cout << "Pruner complete. Likelihood of size 1 at root: " << likelihood[1] << endl;
			for (int i = 0; i<likelihood.size(); ++i)
				cout << "Likelihood of size " << i+1 << " at root: " << likelihood[i] << endl;
		}
        /* END: Computing likelihood of user-specified lambda */

        
        /* START: Estimating lambda(s) (-e) */
        if (estimate) {
			srand(10);
            if (input_file_path.empty()) {
                throw runtime_error("In order to estimate the lambda(s) value(s) (-e), you must specify an input file path (gene family data) with -i. Exiting...");
            }

			if (p_lambda_tree != NULL)
			{
			}


			p_tree->init_gene_family_sizes(*p_gene_families);
			lambda_search_params params(p_tree, *p_gene_families, max_family_size, max_root_family_size);
			single_lambda* p_single_lambda = new single_lambda(find_best_lambda(&params));
			cout << "Best lambda match is " << setw(15) << setprecision(14) << p_single_lambda->get_single_lambda() << endl;

			p_lambda = p_single_lambda;

			return 0;
        }
        /* END: Estimating lambda(s) */
        
        /* START: Running simulations (-s) */
        if (simulate) {
            if (!nsims) {
                throw runtime_error("In order to perform simulations (-s), you must specify the number of simulation runs with -n. Exiting...");
            }

            else { cout << endl << "Performing " << nsims << " simulation batches." << endl; }

            if (input_file_path.empty() && rootdist.empty()) {
                throw runtime_error("In order to perform simulations (s), you must either specify an input file from which the root family size is estimated with -i, or specify a root family distribution with -f. Exiting...");
            }

            /* -i is provided, -f is not */
            else if (rootdist.empty() && !input_file_path.empty()) {
                cout << endl << "Simulations will use the equilibrium root family size estimated from data provided with -i:" << input_file_path << endl;
            }

            /* -f is provided (-f has precedence over -i if both are provided) */
            else {
                cout << "Simulations will use the root family distribution specified with -f: " << rootdist << endl;
                p_rootdist_map = read_rootdist(rootdist); // in map form
                vector<int> rootdist_vec;
                rootdist_vec = vectorize_map(p_rootdist_map); // in vector form

                int max = (*max_element(p_rootdist_map->begin(), p_rootdist_map->end(), max_key<int, int>)).first * 2;

                //cout << "max_family_size = " << max_family_size << endl;
                //vector<trial *> simulation = simulate_families_from_root_size(p_tree, nsims, root_family_size, max_family_size, lambda);
                //vector<vector<trial *> > simulation = simulate_families_from_distribution(p_tree, nsims, *p_rootdist_map, max_family_size, lambda);

                //print_simulation(simulation, cout);

                // total_n_families, lambda_multipliers and lambda_bins will not be harcoded in the future
                int total_n_families = 10;
                int n_cat = 5; // number of gamma categories
                double alpha = 0.5;
                std::vector<double> gamma_cat_probs(n_cat), lambda_multipliers(n_cat);

                get_gamma(gamma_cat_probs, lambda_multipliers, alpha); // passing vectors by reference

                std::vector<int> *p_gamma_cats = weighted_cat_draw(total_n_families, gamma_cat_probs);

                //std::vector<double> lambda_multipliers {1.0, 4.0};
                //std::vector<int> lambda_bins {0, 0, 0, 1, 1, 0, 1, 0, 1, 1}; // the number of elements must be the same as the total key values in p_rootdist_map; here I'm hardcoding it to have 10 elements, as example/test_root_dist.txt has a distribution for 10 families

                rootdist_vec.clear(); // if we want to use uniform (comment to use the file provided with -f)

                // core core_model(cout, &lambda, p_tree, max_family_size, total_n_families, rootdist_vec, lambda_multipliers, *p_gamma_cats);
                // core_model.start_sim_processes();
                // core_model.simulate_processes();
                // core_model.print_simulations(cout);
                // core_model.print_parameter_values();
            }
        }
    
        return 0;
    /* END: Running simulations */
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
