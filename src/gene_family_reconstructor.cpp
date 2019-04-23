#include <cmath>
#include <sstream>

#include "gene_family_reconstructor.h"
#include "lambda.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "gene_family.h"
#include "user_data.h"

std::ostream& operator<<(std::ostream& ost, family_size_change fsc)
{
    switch (fsc)
    {
    case Increase:
        ost << "i";
        break;
    case Decrease:
        ost << "d";
        break;
    case Constant:
        ost << "c";
        break;
    }

    return ost;
}

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

void reconstruct_gene_families(const lambda* lambda, const clade *p_tree,
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
    auto node_id = distance(order.begin(), find(order.begin(), order.end(), node));
    ostringstream ost;
    if (node->is_leaf())
        ost << node->get_taxon_name();
    else
    {
        ost << node_id;
    }

    ost << "_" << textwriter(node);
    if (!node->is_root())
        ost << ':' << node->get_branch_length();
    return ost.str();
}


bool parent_compare(int a, int b)
{
    return a < b;
}

bool parent_compare(double a, double b)
{
    return parent_compare(int(std::round(a)), int(std::round(b)));
}

template <typename T>
void compute_increase_decrease_t(clademap<T>& input, clademap<family_size_change>& output)
{
    for (auto &clade_state : input)
    {
        auto p_clade = clade_state.first;
        T size = clade_state.second;
        if (!p_clade->is_root())
        {
            T parent_size = input[p_clade->get_parent()];
            if (parent_compare(size, parent_size))
                output[p_clade] = Decrease;
            else if (parent_compare(parent_size, size))
                output[p_clade] = Increase;
            else
                output[p_clade] = Constant;
        }
    }
}

void compute_increase_decrease(clademap<int>& input, clademap<family_size_change>& output)
{
    compute_increase_decrease_t(input, output);
}

void compute_increase_decrease(clademap<double>& input, clademap<family_size_change>& output)
{
    compute_increase_decrease_t(input, output);
}

std::ostream& operator<<(std::ostream & ost, const increase_decrease& val)
{
    ost << val.gene_family_id << '\t';
    ost << val.pvalue << "\t";
    ost << (val.pvalue < 0.05 ? 'y' : 'n') << "\t";
    ostream_iterator<family_size_change> out_it(ost, "\t");
    copy(val.change.begin(), val.change.end(), out_it);

    if (!val.category_likelihoods.empty())
    {
        ostream_iterator<double> out_it2(ost, "\t");
        copy(val.category_likelihoods.begin(), val.category_likelihoods.end(), out_it2);
    }
    ost << endl;

    return ost;
}

void reconstruction::write_results(std::string model_identifier, std::string output_prefix, const user_data& data, std::vector<double>& pvalues)
{
    auto order = data.p_tree->find_internal_nodes();

    std::ofstream ofst(filename(model_identifier + "_asr", output_prefix));
    print_reconstructed_states(ofst, order, data.gene_families, data.p_tree);

    std::ofstream family_results(filename(model_identifier + "_family_results", output_prefix));
    print_increases_decreases_by_family(family_results, order, pvalues);

    std::ofstream clade_results(filename(model_identifier + "_clade_results", output_prefix));
    print_increases_decreases_by_clade(clade_results, order);
}

