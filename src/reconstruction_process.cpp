#include <cmath>

#include "reconstruction_process.h"
#include "lambda.h"
#include "io.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"

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


reconstruction_process::reconstruction_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree,
    int max_family_size,
    int max_root_family_size, std::vector<int> rootdist,
    gene_family *gf,
    matrix_cache *p_calc,
    root_equilibrium_distribution* p_prior) : _gene_family(gf), _p_calc(p_calc), _p_prior(p_prior),
    process(ost, lambda, lambda_multiplier, p_tree, max_family_size, max_root_family_size, rootdist)
{
}

class child_multiplier
{
    clademap<std::vector<double> >& _L;
    double value;
    int _j;
public:
    bool write = false;
    child_multiplier(clademap<std::vector<double> >& L, int j) : _L(L), _j(j)
    {
        value = 1.0;
    }

    void operator()(const clade *c)
    {
        value *= _L[c][_j];
        if (write)
            cout << "Child factor " << c->get_taxon_name() << " = " << _L[c][_j] << endl;
    }

    double result() const {
        return value;
    }
};

void reconstruction_process::reconstruct_leaf_node(const clade * c, lambda * _lambda)
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
        L[i] = matrix.get(i, observed_count);
    }
}

void reconstruction_process::reconstruct_root_node(const clade * c)
{
    auto& L = all_node_Ls[c];
    auto& C = all_node_Cs[c];

    L.resize(_max_root_family_size + 1);
    // At the root, we pick a single reconstructed state (step 4 of Pupko)
    C.resize(1);

    // i is the parent, j is the child
    for (size_t i = 1; i < L.size(); ++i)
    {
        double max_val = -1;

        for (size_t j = 1; j < L.size(); ++j)
        {
            child_multiplier cr(all_node_Ls, j);
            c->apply_to_descendants(cr);
            double val = cr.result() * _p_prior->compute(j);
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

void reconstruction_process::reconstruct_internal_node(const clade * c, lambda * _lambda)
{
    auto& C = all_node_Cs[c];
    auto& L = all_node_Ls[c];
    C.resize(_max_family_size + 1);
    L.resize(_max_family_size + 1);

    double branch_length = c->get_branch_length();

    L.resize(_max_family_size + 1);

    auto matrix = _p_calc->get_matrix(branch_length, _lambda->get_value_for_clade(c));

    if (matrix.is_zero())
        throw runtime_error("Zero matrix found");
    // i is the parent, j is the child
    for (size_t i = 0; i < L.size(); ++i)
    {
        size_t max_j;
        double max_val = -1;
        for (size_t j = 0; j < L.size(); ++j)
        {
            child_multiplier cr(all_node_Ls, j);
            c->apply_to_descendants(cr);
            double val = cr.result() * matrix.get(i,j);
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


void reconstruction_process::operator()(const clade *c)
{
    // all_node_Cs and all_node_Ls are hashtables where the keys are nodes and values are vectors of doubles
    unique_ptr<lambda> ml(_lambda->multiply(_lambda_multiplier));

    if (c->is_leaf())
    {
        reconstruct_leaf_node(c, ml.get());
    }
    else if (c->is_root())
    {
        reconstruct_root_node(c);
    }
    else
    {
        reconstruct_internal_node(c, ml.get());
    }
}

class backtracker
{
    std::map<const clade *, std::vector<int> >& _all_node_Cs;
    std::map<const clade *, int>& reconstructed_states;
public:
    backtracker(std::map<const clade *, int>& rc, std::map<const clade *, std::vector<int> >& all_node_Cs, clade *root) : 
        _all_node_Cs(all_node_Cs),
        reconstructed_states(rc)
    {
        reconstructed_states[root] = _all_node_Cs[root][0];
    }

    void operator()(clade *child)
    {
        if (!child->is_leaf())
        {
            auto& C = _all_node_Cs[child];
            int parent_c = reconstructed_states[child->get_parent()];
            reconstructed_states[child] = C[parent_c];
            child->apply_to_descendants(*this);
        }
    }
};

void reconstruction_process::reconstruct()
{
    // Pupko's joint reconstruction algorithm
    _p_tree->apply_reverse_level_order(*this);

    backtracker b(reconstructed_states, all_node_Cs, _p_tree);
    _p_tree->apply_to_descendants(b);

    compute_increase_decrease(reconstructed_states, increase_decrease_map);
}

std::vector<const clade *> reconstruction_process::get_taxa()
{
    vector<const clade *> result;
    for (auto& c : reconstructed_states)
        result.push_back(c.first);

    return result;
}

void reconstruction_process::print_reconstruction(std::ostream & ost, cladevector& order)
{
    ost << _gene_family->id() << '\t';
    for (auto taxon : order)
        ost << reconstructed_states[taxon] << '\t';

    ost << endl;
}

increase_decrease reconstruction_process::get_increases_decreases(cladevector& order, double pvalue)
{
    increase_decrease result;
    result.change.resize(order.size());
    result.gene_family_id = _gene_family->id();
    result.pvalue = pvalue;

    transform(order.begin(), order.end(), result.change.begin(), [this](const clade *taxon)->family_size_change {
        return increase_decrease_map[taxon];
    });

    return result;
}


clademap<double> reconstruction_process::get_weighted_averages(std::vector<reconstruction_process *> m, const vector<double>& _gamma_cat_probs)
{
    vector<const clade *> nodes;
    for (auto& i : m[0]->reconstructed_states)
        nodes.push_back(i.first);

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<_gamma_cat_probs.size(); ++i)
        {
            val += _gamma_cat_probs[i] * double(m[i]->reconstructed_states[node]);
        }
        result[node] = val;
    }

    return result;
}

std::string reconstruction_process::get_family_id() const
{
    return _gene_family->id();
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

void compute_increase_decrease(clademap<int>& input, std::map<const clade *, family_size_change>& output)
{
    compute_increase_decrease_t(input, output);
}

void compute_increase_decrease(clademap<double>& input, std::map<const clade *, family_size_change>& output)
{
    compute_increase_decrease_t(input, output);
}

std::ostream& operator<<(std::ostream & ost, const increase_decrease& val)
{
    ost << val.gene_family_id << '\t';
    ost << (val.pvalue < 0.05 ? 'y' : 'n') << "\t";
    ostream_iterator<family_size_change> out_it(ost, "\t");
    copy(val.change.begin(), val.change.end(), out_it);

    ost << endl;

    return ost;
}

