#include <map>
#include <random>
#include <algorithm>
#include <fstream>

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
    if (argc == 1)
    {
        my_input_parameters.help = true;
        return my_input_parameters;
    }

    int args; // getopt_long returns int or char
    int prev_arg;

    while (prev_arg = optind, (args = getopt_long(argc, argv, "i:e:o:t:y:n:f:E:R:l:m:k:a:s::g::p::r:xb", longopts, NULL)) != -1) {
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
        case 'E':
            my_input_parameters.optimizer_params.neldermead_expansion = atof(optarg);
            break;
        case 'R':
            my_input_parameters.optimizer_params.neldermead_reflection = atof(optarg);
            break;
        case 'f':
            my_input_parameters.rootdist = optarg;
            break;
        case 'r':
            my_input_parameters.chisquare_compare = optarg;
            break;
        case 'g':
            my_input_parameters.do_log = true;
            break;
        case 'h':
            my_input_parameters.help = true;
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

void show_help()
{
    const char *text = ""
        "Usage: cafexp [options]\n\n"
        "CAFE is a software that provides a statistical foundation for evolutionary inferences about changes in gene family size.\n "
        "The program employs a birth and death process to model gene gain and loss across a user-specified phylogenetic tree,\n "
        "thus accounting for the species phylogenetic history. The distribution of family sizes generated under this model can\n "
        "provide a basis for assessing the significance of the observed family size differences among taxa.\n\n"
        "Options:\n"
        "   --fixed_alpha, -a\t\tAlpha value of the discrete gamma distribution to use in category calculations. If not\n \t\t\t\t  specified, the alpha parameter will be estimated by maximum likelihood.\n"
        "   --error_model, -e\t\tError model file path\n"
        "   --rootdist, -f\t\tRoot distribution file path\n"
        "   --infile, -i\t\t\tCharacter or gene family file path\n"
        "   --n_gamma_cats, -k\t\tNumber of gamma rate categories to use. If specified, the Gamma model will be used to run\n \t\t\t\t  calculations, otherwise the Base model will be used.\n"
        "   --fixed_lambda, -l\t\tValue (between 0 and 1) for a single user provided lambda value, otherwise lambda is estimated.\n"
        "   --fixed_multiple_lambdas, -m\tMultiple lambda values, comma separated\n"
        "   --output_prefix, -o\t\tOutput directory - Name of directory automatically created for output\n"
        "   --poisson, -p\t\tUse a Poisson distribution for the root frequency distribution. Without specifying this, a\n \t\t\t\t  normal distribution will be used. A value can be specified (-p10, or --poisson = 10),\n \t\t\t\t  otherwise the distribution will be estimated from the gene families.\n"
        "   --chisquare_compare, -r\tChi square compare\n"
        "   --simulate, -s\t\tSimulate families. Optionally provide the number of simulations to generate (-s100, or --simulate = 100)\n"
        "   --tree, -t\t\t\tTree file path - Required for estimation\n"
        "   --lambda_tree, -y\t\tLambda tree file path\n"
        "   --filter, -x\t\t\tFilter out gene families that don't exist at the root.\n"
        "   --Expansion, -E\t\tExpansion parameter for Nelder-Mead optimizer.\n"
        "   --Reflection, -R\t\tReflection parameter for Nelder-Mead optimizer.\n"
        "   --log, -g\t\t\tName of log file to output.\n"
        "   --lambda_per_family, -b\tEstimate lambda by family (for testing purposes only).\n";

        std::cout << text;
}

/// The main function. Evaluates arguments, calls processes
/// \callgraph
int cafexp(int argc, char *const argv[]) {
    init_lgamma_cache();

    try {
        input_parameters user_input = read_arguments(argc, argv);

        if (user_input.help)
        {
            show_help();
            return 0;
        }
        user_data data;
        data.read_datafiles(user_input);

        if (user_input.exclude_zero_root_families)
        {
            auto rem = std::remove_if(data.gene_families.begin(), data.gene_families.end(), [&data](const gene_family& fam) {
                return !fam.exists_at_root(data.p_tree);
            });

            cout << "\nFiltering the number of families from: " << data.gene_families.size();
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
