#include <map>
#include <random>

#include <getopt.h>

#include "io.h"
#include "execute.h"
#include "simulator.h"

#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "core.h"

using namespace std;

input_parameters read_arguments(int argc, char *const argv[])
{
    input_parameters my_input_parameters;

    int args; // getopt_long returns int or char
    int prev_arg;

    while (prev_arg = optind, (args = getopt_long(argc, argv, "i:e:o:t:y:n:f:l:m:k:a:s::g::p::r:xb", longopts, NULL)) != -1) {
        // while ((args = getopt_long(argc, argv, "i:t:y:n:f:l:e::s::", longopts, NULL)) != -1) {
        if (optind == prev_arg + 2 && optarg && *optarg == '-') {
            cout << "You specified option " << argv[prev_arg] << " but it requires an argument. Exiting..." << endl;
            exit(EXIT_FAILURE);
            // args = ':';
            // --optind;
        }

        switch (args) {
        case 'b':
            my_input_parameters.lambda_per_family = true;
            break;
        case 'i':
            my_input_parameters.input_file_path = optarg;
            break;
        case 'e':
            my_input_parameters.error_model_file_path = optarg;
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
            my_input_parameters.is_simulating = true;
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
        case 'r':
            my_input_parameters.chisquare_compare = optarg;
            break;
        case 'g':
            my_input_parameters.do_log = true;
            break;
        case 'x':
            my_input_parameters.exclude_zero_root_families = true;
            break;
        case ':':   // missing argument
            fprintf(stderr, "%s: option `-%c' requires an argument",
                argv[0], optopt);
            break;
        default: // '?' is parsed (
            throw std::runtime_error(string("Unrecognized parameter: '") + (char)args + "'");
            
        }
    }

    if (optind < argc)
    {
        throw std::runtime_error(string("Unrecognized parameter: '") + argv[optind] + "'");
    }

    my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

    return my_input_parameters;
}

void init_lgamma_cache();

action* get_executor(input_parameters& user_input, user_data& data)
{
    if (!user_input.chisquare_compare.empty()) {
        return new chisquare_compare(data, user_input);
    }
    else if (user_input.is_simulating) {
        return new simulator(data, user_input);
    }
    else
    {
        return new estimator(data, user_input);
    }

    return NULL;
}

/// The main function. Evaluates arguments, calls processes
/// \callgraph
int cafexp(int argc, char *const argv[]) {
    srand(10);
    map<int, int>* p_rootdist_map = NULL; // for sims

    init_lgamma_cache();

    try {
        input_parameters user_input = read_arguments(argc, argv);

        user_data data;
        data.read_datafiles(user_input);

        if (user_input.exclude_zero_root_families)
        {
            auto rem = std::remove_if(data.gene_families.begin(), data.gene_families.end(), [&data](const gene_family& fam) {
                return !fam.exists_at_root(data.p_tree);
            });

            cout << "Filtering the number of families from: " << data.gene_families.size();
            data.gene_families.erase(rem, data.gene_families.end());
            cout << " ==> " << data.gene_families.size() << endl;

        }

        data.p_prior.reset(root_eq_dist_factory(user_input, &data.gene_families));

        // When computing or simulating, only base or gamma model is used. When estimating, base and gamma model are used (to do: compare base and gamma w/ LRT)
        // Build model takes care of -f
        vector<model *> models = build_models(user_input, data);

        unique_ptr<action> act(get_executor(user_input, data));
        if (act)
        {
            act->execute(models);
        }

        /* -g */
        if (user_input.do_log) {

            string prob_matrix_suffix = "_tr_prob_matrices.txt";
            string prob_matrix_file_name = user_input.output_prefix + prob_matrix_suffix;
            std::ofstream ofst(user_input.output_prefix + prob_matrix_suffix);
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
