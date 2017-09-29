#include <cmath>
#include <iomanip>

#include "utils.h"
#include "io.h"
#include "execute.h"
#include "probability.h"
#include "lambda.h"
#include "core.h"
#include "gamma_core.h"

//! Read user provided gene family data (whose path is stored in input_parameters instance)
std::vector<gene_family> * execute::read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree) {
    
    std::vector<gene_family> *p_gene_families = NULL;
    
    if (!my_input_parameters.input_file_path.empty()) {
        ifstream input_file(my_input_parameters.input_file_path); 

        p_gene_families = read_gene_families(input_file, p_tree);
            
        // Iterating over gene families to get max gene family size
        for (std::vector<gene_family>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
            int this_family_max_size = it->get_max_size();
            
            if (max_family_size < this_family_max_size)
                max_family_size = this_family_max_size;
        }

        max_root_family_size = std::max(30, static_cast<int>(std::rint(max_family_size*1.25)));
        max_family_size = max_family_size + std::max(50, max_family_size/5);
        cout << "Max family size is: " << max_family_size << endl;
        cout << "Max root family size is: " << max_root_family_size << endl;
            
        }

    else { p_gene_families = new std::vector<gene_family>(); }
    
    return p_gene_families;
}

//! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * execute::read_input_tree(const input_parameters &my_input_parameters) {
    clade *p_tree = new clade();
    p_tree = read_tree(my_input_parameters.tree_file_path, false);
    
    return p_tree;
}

//! Read user provided lambda tree (lambda structure)
clade * execute::read_lambda_tree(const input_parameters &my_input_parameters) {
    clade * p_lambda_tree = NULL; // not new clade() because read_tree() below also allocates memory (so this would leak memory)
    
    if (!my_input_parameters.lambda_tree_file_path.empty()) {
        p_lambda_tree = read_tree(my_input_parameters.lambda_tree_file_path, true);
        // node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
        // p_lambda_tree->print_clade();
        }
    
    return p_lambda_tree;
}

//! Read user provided single or multiple lambdas
lambda * execute::read_lambda(const input_parameters &my_input_parameters, probability_calculator &my_calculator, clade * p_lambda_tree) {
      
    lambda *p_lambda = NULL; // lambda is an abstract class, and so we can only instantiate it as single_lambda or multiple lambda -- therefore initializing it to NULL
        
    // -l
    if (my_input_parameters.fixed_lambda > 0.0) {
        cout << "Specified lambda (-l): " << my_input_parameters.fixed_lambda << endl;
        p_lambda = new single_lambda(&my_calculator, my_input_parameters.fixed_lambda);
            // vector<double> posterior = get_posterior((*p_gene_families), max_family_size, fixed_lambda, p_tree);
            // double map = log(*max_element(posterior.begin(), posterior.end()));
            // cout << "Posterior values found - max log posterior is " << map << endl;
			// call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
        }
    
    // -m
    if (!my_input_parameters.fixed_multiple_lambdas.empty()) {       
        if (p_lambda_tree == NULL) {
            throw runtime_error("You must specify a lambda tree (-y) if you fix multiple lambda values (-m). Exiting..."); 
        }
        
        map<std::string, int> node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map(); // allows matching lambda values (when -m) to clades in lambda tree

        vector<string> lambdastrings = tokenize_str(my_input_parameters.fixed_multiple_lambdas, ',');
	vector<double> lambdas(lambdastrings.size());

        // transform is like R's apply (lambdas takes the outputs)
        transform(lambdastrings.begin(), lambdastrings.end(), lambdas.begin(),
                [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
                );
        
        p_lambda = new multiple_lambda(&my_calculator, node_name_to_lambda_index, lambdas);
    }
    
    return p_lambda;
} 

std::string filename(std::string base, std::string suffix)
{
    return base + (suffix.empty() ? "" : "_") + suffix + ".txt";
}

void execute::infer(std::vector<core *>& models, std::vector<gene_family> *p_gene_families, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size)
{
    std::ofstream results(filename("results", my_input_parameters.output_prefix));

    std::ofstream likelihoods(filename("family_lks", my_input_parameters.output_prefix));

    for (int i = 0; i < models.size(); ++i) {
        models[i]->set_gene_families(p_gene_families);
        
        gamma_core* p_model = dynamic_cast<gamma_core *>(models[i]);
        if (p_model != NULL) {
            p_model->initialize_with_alpha(my_input_parameters.n_gamma_cats, p_gene_families->size(), my_input_parameters.fixed_alpha);
        }

        models[i]->set_max_sizes(max_family_size, max_root_family_size);

        cout << endl << "Starting inference processes for model " << i << endl;
        models[i]->start_inference_processes();

        cout << endl << "Inferring processes for model " << i << endl;
        double result = models[i]->infer_processes();
        results << "Model " << models[i]->name() << " Result: " << result << endl;

        likelihoods << "#FamilyID\tGamma Cat Median\tLikelihood of Category\tLikelihood of Family" << endl;
        models[i]->print_results(likelihoods);
    }
}

void execute::simulate(std::vector<core *>& models, const input_parameters &my_input_parameters)
{
    //                cout << "Simulations will use the root family distribution specified with -f: " << my_input_parameters.rootdist << endl;
    std::map<int, int> *p_rootdist_map = read_rootdist(my_input_parameters.rootdist); // in map form
    vector<int> rootdist_vec;

    if (!p_rootdist_map)
        throw runtime_error("p_rootdist_map not specified. Exiting...");

    rootdist_vec = vectorize_map(p_rootdist_map); // in vector form
    rootdist_vec.clear();

    cout << "I will simulate with this many models: " << models.size() << endl;

    for (int i = 0; i < models.size(); ++i) {

        cout << "Simulating for model " << i << endl;

        models[i]->set_rootdist_vec(rootdist_vec);

        cout << "I just set the root distribution." << endl;

        // max_family_size = (*max_element(p_rootdist_map->begin(), p_rootdist_map->end(), max_key<int, int>)).first * 2;

        // cout << "My max family size is: " << max_family_size << " and my max root family size is: " << max_root_family_size << endl;

        //cout << "max_family_size = " << max_family_size << endl;
        //vector<trial *> simulation = simulate_families_from_root_size(p_tree, nsims, root_family_size, max_family_size, lambda);
        //vector<vector<trial *> > simulation = simulate_families_from_distribution(p_tree, nsims, *p_rootdist_map, max_family_size, lambda);

        //print_simulation(simulation, cout);

        // lambda_multipliers and lambda_bins will not be harcoded in the future
        models[i]->set_total_n_families_sim(my_input_parameters.nsims);

        cout << "I just set number of families to simulate." << endl;

        gamma_core* p_model = dynamic_cast<gamma_core *>(models[i]);
        if (p_model != NULL) {
            bool has_alpha = false;
            if (has_alpha)
            {
                p_model->initialize_with_alpha(my_input_parameters.n_gamma_cats, my_input_parameters.nsims, 0.5);
            }
            else
            {
                vector<double> lambda_multipliers = { 1.0, 4.0 };
                std::vector<int> gamma_cats{ 0, 0, 0, 1, 1, 0, 1, 0, 1, 1 }; // the number of elements must be the same as the total key values in p_rootdist_map; here I'm hardcoding it to have 10 elements, as example/test_root_dist.txt has a distribution for 10 families
                p_model->initialize_without_alpha(my_input_parameters.n_gamma_cats, my_input_parameters.nsims,
                    lambda_multipliers, gamma_cats);
            }
        }

        // core core_model(cout, p_lambda, p_tree, max_family_size, total_n_families, rootdist_vec, n_cat, alpha);

        // model(cout, p_lambda, p_tree, max_family_size, total_n_families, rootdist_vec, lambda_bins, lambda_multipliers);
        models[i]->start_sim_processes();
        
        std::ofstream ofst(filename("simulation", my_input_parameters.output_prefix));
        models[i]->simulate_processes();
        models[i]->print_processes(ofst);
        //    for (trial::iterator it = _my_simulation->begin(); it != _my_simulation->end(); ++it) {
        //        ost << "#" << it->first->get_taxon_name() << "\n";
        //    }


        // models[i]->adjust_family(cout);
        // core_model.print_parameter_values();
    }
}

lambda * execute::estimate_lambda(const input_parameters &my_input_parameters, clade *p_tree, clade *p_lambda_tree, 
    std::vector<gene_family>* p_gene_families, int max_family_size, int max_root_family_size, probability_calculator& calculator)
{
    if (my_input_parameters.input_file_path.empty()) {
        throw runtime_error("In order to estimate the lambda(s) value(s) (-e), you must specify an input file path (gene family data) with -i. Exiting...");
    }

    if (p_lambda_tree != NULL)
    {
    }


    p_tree->init_gene_family_sizes(*p_gene_families);
    lambda_search_params params(p_tree, *p_gene_families, max_family_size, max_root_family_size);
    single_lambda* p_single_lambda = new single_lambda(&calculator, find_best_lambda(&params));
    cout << "Best lambda match is " << setw(15) << setprecision(14) << p_single_lambda->get_single_lambda() << endl;

    return p_single_lambda;
}