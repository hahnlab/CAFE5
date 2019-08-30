#include <cmath>
#include <sstream>

#include "gene_family_reconstructor.h"
#include "lambda.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "gene_family.h"
#include "user_data.h"

void reconstruct_leaf_node(const clade * c, const lambda * _lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int _max_family_size, const gene_family* _gene_family, const matrix_cache *_p_calc)
{
    auto& C = all_node_Cs[c];
    auto& L = all_node_Ls[c];
    C.resize(_max_family_size + 1);
    L.resize(_max_family_size + 1);

    double branch_length = c->get_branch_length();

    L.resize(_max_family_size + 1);

    int observed_count = _gene_family->get_species_size(c->get_taxon_name());
    fill(C.begin(), C.end(), observed_count);

    auto matrix = _p_calc->get_matrix(branch_length, _lambda->get_value_for_clade(c));
    // i will be the parent size
    for (size_t i = 1; i < L.size(); ++i)
    {
        L[i] = matrix->get(i, observed_count);
    }
}

void reconstruct_root_node(const clade * c, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int _max_family_size, int _max_root_family_size, const root_equilibrium_distribution* _p_prior)
{
    auto& L = all_node_Ls[c];
    auto& C = all_node_Cs[c];

    L.resize(min(_max_family_size, _max_root_family_size) + 1);
    // At the root, we pick a single reconstructed state (step 4 of Pupko)
    C.resize(1);

    // i is the parent, j is the child
    for (size_t i = 1; i < L.size(); ++i)
    {
        double max_val = -1;

        for (size_t j = 1; j < L.size(); ++j)
        {
            double value = 1.0;
            auto child_multiplier = [&all_node_Ls, j, &value](const clade *child) {
                value *= all_node_Ls[child][j];
            };
            c->apply_to_descendants(child_multiplier);
            double val = value * _p_prior->compute(j);
            if (val > max_val)
            {
                max_val = val;
                C[0] = j;
            }
        }

        L[i] = max_val;
    }

    
    if (*max_element(L.begin(), L.end()) == 0.0)
    {
        cerr << "WARNING: failed to calculate L value at root" << endl;
    }
}

void reconstruct_internal_node(const clade * c, const lambda * _lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int _max_family_size, const matrix_cache *_p_calc)
{
    auto& C = all_node_Cs[c];
    auto& L = all_node_Ls[c];
    C.resize(_max_family_size + 1);
    L.resize(_max_family_size + 1);

    double branch_length = c->get_branch_length();

    L.resize(_max_family_size + 1);

    auto matrix = _p_calc->get_matrix(branch_length, _lambda->get_value_for_clade(c));

    if (matrix->is_zero())
        throw runtime_error("Zero matrix found");
    // i is the parent, j is the child
    for (size_t i = 0; i < L.size(); ++i)
    {
        size_t max_j = 0;
        double max_val = -1;
        for (size_t j = 0; j < L.size(); ++j)
        {
            double value = 1.0;
            auto child_multiplier = [&all_node_Ls, j, &value](const clade *child) {
                value *= all_node_Ls[child][j];
            };
            c->apply_to_descendants(child_multiplier);
            double val = value * matrix->get(i,j);
            if (val > max_val)
            {
                max_j = j;
                max_val = val;
            }
        }

        L[i] = max_val;
        C[i] = max_j;
    }
}


void reconstruct_at_node(const clade *c, const lambda *_lambda, clademap<std::vector<int>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int max_family_size, int max_root_family_size, const matrix_cache* p_calc, const root_equilibrium_distribution* p_prior, const gene_family *p_family)
{
    if (c->is_leaf())
    {
        reconstruct_leaf_node(c, _lambda, all_node_Cs, all_node_Ls, max_family_size, p_family, p_calc);
    }
    else if (c->is_root())
    {
        reconstruct_root_node(c, all_node_Cs, all_node_Ls, max_family_size, max_root_family_size, p_prior);
    }
    else
    {
        reconstruct_internal_node(c, _lambda, all_node_Cs, all_node_Ls, max_family_size, p_calc);
    }
}

void reconstruct_gene_family(const lambda* lambda, const clade *p_tree,
    int max_family_size,
    int max_root_family_size,
    const gene_family *gf,
    matrix_cache *p_calc,
    root_equilibrium_distribution* p_prior, clademap<int>& reconstructed_states)
{
    clademap<std::vector<int>> all_node_Cs;

    /// Ls hold a probability for each family size (values are probabilities of any given family size)
    clademap<std::vector<double>> all_node_Ls;

    std::function <void(const clade *)> pupko_reconstructor;
    pupko_reconstructor = [&](const clade *c) {
        reconstruct_at_node(c, lambda, all_node_Cs, all_node_Ls, max_family_size, max_root_family_size, p_calc, p_prior, gf);
    };

    std::function<void(const clade *child)> backtracker;
    backtracker = [&reconstructed_states, &all_node_Cs, &backtracker](const clade *child) {
        if (!child->is_leaf())
        {
            auto& C = all_node_Cs[child];
            int parent_c = reconstructed_states[child->get_parent()];
            reconstructed_states[child] = C[parent_c];
            child->apply_to_descendants(backtracker);
        }
        };

    // Pupko's joint reconstruction algorithm
    p_tree->apply_reverse_level_order(pupko_reconstructor);

    reconstructed_states[p_tree] = all_node_Cs[p_tree][0];
    p_tree->apply_to_descendants(backtracker);

}

string newick_node(const clade *node, const cladevector& order, std::function<std::string(const clade *c)> textwriter)
{
    ostringstream ost;
    ost << clade_index_or_name(node, order) << "_" << textwriter(node);

    if (!node->is_root())
        ost << ':' << node->get_branch_length();

    return ost.str();
}


int parent_compare(int a, int b)
{
    return a - b;
}

int parent_compare(double a, double b)
{
    return parent_compare(int(std::round(a)), int(std::round(b)));
}

template <typename T>
void compute_increase_decrease_t(clademap<T>& input, clademap<int>& output)
{
    for (auto &clade_state : input)
    {
        auto p_clade = clade_state.first;
        T size = clade_state.second;
        if (!p_clade->is_root())
        {
            T parent_size = input[p_clade->get_parent()];
            output[p_clade] = parent_compare(size, parent_size);
        }
        else
        {
            output[p_clade] = 0;
        }
    }
}

void compute_increase_decrease(clademap<int>& input, clademap<int>& output)
{
    compute_increase_decrease_t(input, output);
}

void compute_increase_decrease(clademap<double>& input, clademap<int>& output)
{
    compute_increase_decrease_t(input, output);
}

void reconstruction::print_increases_decreases_by_family(std::ostream& ost, const cladevector& order, const std::vector<const gene_family*>& gene_families, const std::vector<double>& pvalues) {
    if (gene_families.size() != pvalues.size())
    {
        throw std::runtime_error("No pvalues found for family");
    }
    if (gene_families.empty())
    {
        ost << "No increases or decreases recorded\n";
        return;
    }

    ost << "#FamilyID\tpvalue\t*\t";
    for (auto& it : order) {
        ost << clade_index_or_name(it, order) << "\t";
    }
    ost << endl;

    for (size_t i = 0; i < gene_families.size(); ++i) {
        ost << gene_families[i]->id() << '\t' << pvalues[i] << '\t';
        ost << (pvalues[i] < 0.05 ? 'y' : 'n');
        for (auto c : order)
            ost << '\t' << get_increase_decrease(gene_families[i], c) ;
        ost << endl;
    }
}

void reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order, const std::vector<const gene_family*>& gene_families) {
    clademap<pair<int, int>> increase_decrease_map;

    for (size_t j = 0; j < gene_families.size(); ++j) {
        for (size_t i = 0; i < order.size(); ++i)
        {
            int val = get_delta(gene_families[j], order[i]);
            if (val > 0)
                increase_decrease_map[order[i]].first++;
            if (val < 0)
                increase_decrease_map[order[i]].second++;
        }
    }

    ost << "#Taxon_ID\t#Taxon_Name\tIncrease/Decrease\n";
    for (auto& it : increase_decrease_map) {
        ost << clade_index_or_name(it.first, order) << "\t";
        ost << it.first->get_taxon_name() << "\t";
        ost << it.second.first << "/" << it.second.second << endl;
    }
}

void reconstruction::print_family_clade_table(std::ostream& ost, const cladevector& order, const std::vector<const gene_family*>& gene_families, const clade* p_tree, std::function<string(int family_index, const clade *c)> get_family_clade_value)
{
    ost << "FamilyID";
    for (auto c : order)
    {
        ost << "\t" << clade_index_or_name(c, order);
    }
    ost << endl;
    for (size_t i = 0; i < gene_families.size(); ++i)
    {
        ost << gene_families[i]->id();
        for (auto node : order)
        {
            ost << "\t";
            ost << get_family_clade_value(i, node);
        }
        ost << endl;
    }
}

void print_branch_probabilities(std::ostream& ost, const cladevector& order, const std::vector<const gene_family*>& gene_families, const vector<clademap<double>>& branch_probabilities)
{
    ost << "#FamilyID\t";
    for (auto& it : order) {
        ost << clade_index_or_name(it, order) << "\t";
    }
    ost << endl;

    for (size_t i = 0; i < gene_families.size(); ++i) {
        ost << gene_families[i]->id();
        for (auto c : order)
            ost << '\t' << branch_probabilities[i].at(c);
        ost << endl;
    }

}

void reconstruction::write_results(std::string model_identifier, std::string output_prefix, const clade *p_tree, const std::vector<const gene_family*>& families, std::vector<double>& pvalues, const std::vector<clademap<double>>& branch_probabilities)
{
    cladevector order;
    auto fn = [&order](const clade* c) { order.push_back(c);  };
    p_tree->apply_reverse_level_order(fn);

    std::ofstream ofst(filename(model_identifier + "_asr", output_prefix));
    print_reconstructed_states(ofst, order, families, p_tree);

    std::ofstream counts(filename(model_identifier, output_prefix) + ".count");
    print_node_counts(counts, order, families, p_tree);

    std::ofstream change(filename(model_identifier, output_prefix) + ".change");
    print_node_change(change, order, families, p_tree);

    std::ofstream family_results(filename(model_identifier + "_family_results", output_prefix));
    print_increases_decreases_by_family(family_results, order, families, pvalues);

    std::ofstream clade_results(filename(model_identifier + "_clade_results", output_prefix));
    print_increases_decreases_by_clade(clade_results, order, families);

    std::ofstream branch_probabilities_file(filename(model_identifier + "_branch_probabilities", output_prefix));
    print_branch_probabilities(branch_probabilities_file, order, families, branch_probabilities);

    print_additional_data(order, families, output_prefix);
}

void viterbi_sum_probabilities(const clade* c, const gene_family& family, const reconstruction* rec, int max_family_size, const matrix_cache& cache, const lambda* p_lambda, clademap<double>& results)
{
    if (c->is_root())
    {
        results[c] = 0;
        return;
    }

    const matrix* probs = cache.get_matrix(c->get_branch_length(), p_lambda->get_value_for_clade(c));

    int parent_size = rec->reconstructed_size(family, c->get_parent());
    double calculated_probability = probs->get(parent_size, rec->reconstructed_size(family, c));
    for (int m = 0; m < max_family_size; m++)
    {
        double probability_to_m = probs->get(parent_size, m);
        if (probability_to_m == calculated_probability)
        {
            results[c] += probability_to_m / 2.0;
        }
        else if (probability_to_m < calculated_probability)
        {
            results[c] += probability_to_m;
        }
    }
}

clademap<double> compute_branch_level_probabilities(const clade* p_tree, const gene_family& family, const reconstruction* rec, const lambda* p_lambda, const matrix_cache& cache, int max_family_size, int max_root_family_size)
{
    clademap<double> results;

    auto sum_probabilities_func = [&](const clade* parent) {
        viterbi_sum_probabilities(parent, family, rec, max_family_size, cache, p_lambda, results);
    };
    p_tree->apply_reverse_level_order(sum_probabilities_func);

    return results;
}
