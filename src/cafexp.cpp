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
#include "chisquare.h"
#include "matrix_cache.h"

using namespace std;

double pvalue(double v, const vector<double>& conddist)
{
    int idx = conddist.size() - 1;

    auto bound = std::upper_bound(conddist.begin(), conddist.end(), v);
    if (bound != conddist.end())
    {
        idx = bound - conddist.begin();
    }
    return  idx / (double)conddist.size();
}


// find_fast_families under base model through simulations (if we reject gamma)
vector<double> compute_pvalues(int max_family_size, int max_root_family_size, int number_of_simulations, lambda *p_lambda, vector<gene_family>& families, clade* p_tree)
{
    cout << "Computing pvalues\n" << std::unitbuf;

	double lambda_val = dynamic_cast<single_lambda *>(p_lambda)->get_single_lambda();
    matrix_cache cache;
    branch_length_finder lengths;
    p_tree->apply_prefix_order(lengths);
    cache.precalculate_matrices(max_family_size + 1, get_lambda_values(p_lambda), lengths.result());
    
    auto cd = get_conditional_distribution_matrix(p_tree, max_root_family_size, max_family_size, number_of_simulations, lambda_val, cache);

    vector<double> result(families.size());
    transform(families.begin(), families.end(), result.begin(), [max_family_size, p_lambda, p_tree, &cache, &cd](gene_family& gf)->double {
        cout << ".";
		int max = gf.get_max_size();
		double max_root_family_size = rint(max*1.25);

        map<string, int> species_count;
        for (auto& species : gf.get_species()) {
            species_count[species] = gf.get_species_size(species);
        }

		likelihood_computer pruner(max_root_family_size, max_family_size, p_lambda, species_count, cache, NULL);
		p_tree->apply_reverse_level_order(pruner);
		auto lh = pruner.get_likelihoods(p_tree);

		double observed_max_likelihood = pruner.max_likelihood(p_tree);	// max value but do we need a posteriori value instead?

        vector<double> pvalues(max_root_family_size);
		for (int s = 0; s < max_root_family_size; s++)
		{
			pvalues[s] = pvalue(observed_max_likelihood, cd.at(s));
		}
        return *max_element(pvalues.begin(), pvalues.end());
    });

    cout << "done!\n";

    return result;
}

input_parameters read_arguments(int argc, char *const argv[])
{
    input_parameters my_input_parameters;

    int args; // getopt_long returns int or char
    int prev_arg;

    while (prev_arg = optind, (args = getopt_long(argc, argv, "i:e:o:t:y:n:f:l:m:k:a:s::g::p::r:", longopts, NULL)) != -1) {
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
        case ':':   // missing argument
            fprintf(stderr, "%s: option `-%c' requires an argument",
                argv[0], optopt);
            break;
        default: // '?' is parsed (
            throw std::runtime_error(string("Unrecognized parameter: '") + (char)args + "'");
            
        }
    }

    my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

    return my_input_parameters;
}

void show_pvalues(string values)
{
    vector<string> chistrings = tokenize_str(values, ',');
    vector<double> chis(chistrings.size());

    // transform is like R's apply (vector lambdas takes the outputs, here we are making doubles from strings
    transform(chistrings.begin(), chistrings.end(), chis.begin(),
        [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
    );

    double degrees_of_freedom = chis[2];
    cout << "PValue = " << 1.0 - chi2cdf(2 * (chis[1] - chis[0]), degrees_of_freedom) << std::endl;
}

struct user_data {
    int max_family_size = -1; // needed for defining matrix size
    int max_root_family_size = -1; // needed for defining matrix size

    clade *p_tree = NULL; // instead of new clade(), o.w. mem leak
    lambda *p_lambda;
    clade *p_lambda_tree = NULL;
    error_model *p_error_model = NULL;
    std::vector<gene_family> gene_families;
};

user_data read_datafiles(execute& my_executer, const input_parameters& my_input_parameters)
{
    user_data data;

    /* -t */
    if (!my_input_parameters.tree_file_path.empty()) {
        data.p_tree = my_executer.read_input_tree(my_input_parameters); // populates p_tree (pointer to phylogenetic tree)
    }

    /* -i */
    if (!my_input_parameters.input_file_path.empty()) {
        // Populates (pointer to) vector of gene family instances, max_family_size and max_root_family_size (last two passed by reference)
        my_executer.read_gene_family_data(my_input_parameters, data.max_family_size, data.max_root_family_size, data.p_tree, &data.gene_families);
    }

    /* -e */
    if (!my_input_parameters.error_model_file_path.empty()) {
        data.p_error_model = new error_model;
        my_executer.read_error_model(my_input_parameters, data.p_error_model);
    }

    /* -y */
    if (!my_input_parameters.lambda_tree_file_path.empty()) {
        data.p_lambda_tree = my_executer.read_lambda_tree(my_input_parameters);
    }

    /* -l/-m (in the absence of -l, estimate) */
    data.p_lambda = my_executer.read_lambda(my_input_parameters, data.p_lambda_tree);

    return data;
}

/// The main function. Evaluates arguments, calls processes
/// \callgraph
int cafexp(int argc, char *const argv[]) {
    /* START: Option variables for main() */
    srand(10);
    map<int, int>* p_rootdist_map = NULL; // for sims

    execute my_executer;
    try {
        input_parameters my_input_parameters = read_arguments(argc, argv);

        /* -r */
        if (!my_input_parameters.chisquare_compare.empty()) {

            show_pvalues(my_input_parameters.chisquare_compare);

            return 0;
        }

        user_data data = read_datafiles(my_executer, my_input_parameters);

        root_equilibrium_distribution* p_prior = root_eq_dist_factory(my_input_parameters, &data.gene_families);
        
        // When computing or simulating, only base or gamma model is used. When estimating, base and gamma model are used (to do: compare base and gamma w/ LRT)
        // Build model takes care of -f
        vector<model *> models = build_models(my_input_parameters, data.p_tree, data.p_lambda, &data.gene_families, data.max_family_size, data.max_root_family_size, data.p_error_model);
        
        if (models.empty())
            throw std::runtime_error("Not enough information to specify a model");
        
        /* -s */
        if (my_input_parameters.is_simulating()) {
            // -s is provided an argument (-f is not), using -i to obtain root eq freq distr'n
            if (my_input_parameters.nsims != 0) {
                throw std::runtime_error("A specified number of simulations with a root distribution is not supported");
                // place holder for estimating poisson lambda if -p, or using uniform as root eq freq distr'n
            }

            // -f is provided (-s does not have an argument), not using -i
            else if (!my_input_parameters.rootdist.empty()) {
                cout << "Using -f, not using -i, nsims = " << my_input_parameters.nsims << endl;
                my_executer.simulate(models, my_input_parameters);
            }
        }
        else {  // not simulating - calculate values
            if (data.p_lambda == NULL)
            {
                my_executer.estimate_lambda(models, data.gene_families, data.p_error_model, data.p_tree, data.p_lambda_tree, p_prior);
                data.p_lambda = models[0]->get_lambda();
            }

            my_executer.compute(models, &data.gene_families, p_prior, my_input_parameters, data.max_family_size, data.max_root_family_size);

            my_executer.reconstruct(models, my_input_parameters, p_prior);

            auto pvalues = compute_pvalues(data.max_family_size, data.max_root_family_size, 1000, data.p_lambda, data.gene_families, data.p_tree);

            my_executer.write_results(models, my_input_parameters, pvalues);

        }
        delete p_prior;


        /* -g */
        if (my_input_parameters.do_log) {

            string prob_matrix_suffix = "_tr_prob_matrices.txt";
            string prob_matrix_file_name = my_input_parameters.output_prefix + prob_matrix_suffix;
            std::ofstream ofst(my_input_parameters.output_prefix + prob_matrix_suffix);
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
