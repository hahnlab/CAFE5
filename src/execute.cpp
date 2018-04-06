#include <cmath>
#include <set>

#include "utils.h"
#include "io.h"
#include "execute.h"
#include "probability.h"
#include "lambda.h"
#include "core.h"
#include "gamma_core.h"
#include "root_equilibrium_distribution.h"
#include "chisquare.h"
#include "matrix_cache.h"
#include "base_model.h"

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


//! Read user provided gene family data (whose path is stored in input_parameters instance)
void execute::read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree, std::vector<gene_family> *p_gene_families) {
    
    ifstream input_file(my_input_parameters.input_file_path); 
    if (!input_file.is_open())
        throw std::runtime_error("Failed to open " + my_input_parameters.input_file_path + ". Exiting...");
    
    read_gene_families(input_file, p_tree, p_gene_families); // in io.cpp/io.h
            
    // Iterating over gene families to get max gene family size
    for (std::vector<gene_family>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
        int this_family_max_size = it->get_max_size();
            
        if (max_family_size < this_family_max_size)
            max_family_size = this_family_max_size;
    }

    max_root_family_size = std::max(30, static_cast<int>(std::rint(max_family_size*1.25)));
    max_family_size = max_family_size + std::max(50, max_family_size/5);
    // cout << "Read input file " << my_input_parameters.input_file_path << "." << endl;
    // cout << "Max (parsed) family size is: " << max_family_size << endl;
    // cout << "Max root family size is: " << max_root_family_size << endl;
}

//! Read user provided error model file (whose path is stored in input_parameters instance)
void execute::read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model) {
    
    ifstream error_model_file(my_input_parameters.error_model_file_path); 
    if (!error_model_file.is_open()) {
        throw std::runtime_error("Failed to open " + my_input_parameters.error_model_file_path + ". Exiting...");
    }

    read_error_model_file(error_model_file, p_error_model);
    
} // GOTTA WRITE THIS!

//! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * execute::read_input_tree(const input_parameters &my_input_parameters) {    
    return read_tree(my_input_parameters.tree_file_path, false);
}

//! Read user provided lambda tree (lambda structure)
clade * execute::read_lambda_tree(const input_parameters &my_input_parameters) {
        return read_tree(my_input_parameters.lambda_tree_file_path, true);
}

//! Read user provided single or multiple lambdas
lambda * execute::read_lambda(const input_parameters &my_input_parameters, clade *p_lambda_tree) {
      
    lambda *p_lambda = NULL; // lambda is an abstract class, and so we can only instantiate it as single_lambda or multiple lambda -- therefore initializing it to NULL
        
    // -l
    if (my_input_parameters.fixed_lambda > 0.0) {
        p_lambda = new single_lambda(my_input_parameters.fixed_lambda);
        // call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
    }
    
    // -m
    if (!my_input_parameters.fixed_multiple_lambdas.empty()) {       
        map<std::string, int> node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map(); // allows matching different lambda values to nodes in lambda tree
        vector<string> lambdastrings = tokenize_str(my_input_parameters.fixed_multiple_lambdas, ',');
	vector<double> lambdas(lambdastrings.size());

        // transform is like R's apply (vector lambdas takes the outputs, here we are making doubles from strings
        transform(lambdastrings.begin(), lambdastrings.end(), lambdas.begin(),
                [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
                );
        
        p_lambda = new multiple_lambda(node_name_to_lambda_index, lambdas);
    }
    
    return p_lambda;
} 

std::string filename(std::string base, std::string suffix)
{
    return base + (suffix.empty() ? "" : "_") + suffix + ".txt";
}

void execute::compute(std::vector<model *>& models, std::vector<gene_family> *p_gene_families, root_equilibrium_distribution *p_prior, const input_parameters &my_input_parameters, int max_family_size, int max_root_family_size)
{
    std::ofstream results_file(filename("results", my_input_parameters.output_prefix));

    std::ofstream likelihoods_file(filename("family_lks", my_input_parameters.output_prefix));

    std::vector<double> model_likelihoods(models.size());
    for (int i = 0; i < models.size(); ++i) {
        cout << endl << "Starting inference processes for " << models[i]->name() << " model" << endl;
        models[i]->start_inference_processes();

        cout << endl << "Inferring processes for " << models[i]->name() << " model" << endl;
        double result = models[i]->infer_processes(p_prior);
        models[i]->write_vital_statistics(results_file, result);

        models[i]->write_family_likelihoods(likelihoods_file);

        model_likelihoods[i] = result;
    }

    if (model_likelihoods.size() == 2)
    {
        cout << "PValue = " << (1.0 - chi2cdf(2 * (model_likelihoods[1] - model_likelihoods[0]), 1.0));
    }
}

bool compare_result(const optimizer::result& a, const optimizer::result& b)
{
    return a.score < b.score;
}

void execute::estimate_lambda(std::vector<model *>& models, std::vector<gene_family> &gene_families, error_model *p_error_model, clade *p_tree, clade *p_lambda_tree, root_equilibrium_distribution *p_prior)
{
    if (p_tree == NULL)
    {
        throw runtime_error("No tree specified for lambda estimation");
    }
    if (p_error_model)
    {   // we only support base model epsilon optimizing at the moment
        base_model *b = dynamic_cast<base_model *>(models[0]);
        b->initialize_lambda(p_lambda_tree);
        unique_ptr<optimizer_scorer> scorer(b->get_epsilon_optimizer(p_prior));
        optimizer opt(scorer.get());
        vector<optimizer::result> results;
        for (double epsilon = 0.05; epsilon < .5; epsilon += .1)
        {
            p_error_model->update_single_epsilon(epsilon);
            results.push_back(opt.optimize());
        }
        auto best = min_element(results.begin(), results.end(), compare_result);
        scorer->finalize(&best->values[0]);
        cout << "Final score: " << best->score << ", Lambda: " << best->values[0] << ", Epsilon: " << best->values[1] * 2 << endl;
    }
    else
    {
        for (model* p_model : models) {
            p_model->initialize_lambda(p_lambda_tree);

            p_tree->init_gene_family_sizes(gene_families);
            unique_ptr<optimizer_scorer> scorer(p_model->get_lambda_optimizer(p_prior));
            optimizer opt(scorer.get());
            auto result = opt.optimize();
            scorer->finalize(&result.values[0]);
        }
    }

}

/// Simulate
/// \callgraph
void execute::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    cout << "I will simulate with this many models: " << models.size() << endl;

    for (int i = 0; i < models.size(); ++i) {

        cout << "Simulating for model " << models[i]->name() << endl;

        // lambda_multipliers and lambda_bins will not be harcoded in the future
        if (my_input_parameters.rootdist.empty())
            models[i]->set_total_n_families_sim(my_input_parameters.nsims);

        cerr << "Simulating " << models[i]->get_total_n_families_sim() << " families" << endl;

        models[i]->start_sim_processes();
        
        std::ofstream ofst(filename("simulation", my_input_parameters.output_prefix));
        models[i]->simulate_processes();
        models[i]->print_processes(ofst);
    }
}

void execute::reconstruct(std::vector<model *>& models, const input_parameters &my_input_parameters, root_equilibrium_distribution *p_prior)
{
    matrix_cache cache;
    for (model* p_model : models) {
        p_model->reconstruct_ancestral_states(&cache, p_prior);
    }
}

void execute::write_results(std::vector<model *>& models, const input_parameters &my_input_parameters, std::vector<double>& pvalues)
{
    for (model* p_model : models) {
        std::ofstream ofst(filename(p_model->name() + "_asr", my_input_parameters.output_prefix));
        p_model->print_reconstructed_states(ofst);

        std::ofstream family_results(filename(p_model->name() + "_family_results", my_input_parameters.output_prefix));
        p_model->print_increases_decreases_by_family(family_results, pvalues);

        std::ofstream clade_results(filename(p_model->name() + "_clade_results", my_input_parameters.output_prefix));
        p_model->print_increases_decreases_by_clade(clade_results);
    }

    std::ofstream dfst(filename("depth", my_input_parameters.output_prefix));
    models[0]->print_node_depths(dfst);
}
