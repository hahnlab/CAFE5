#include "reconstruction_process.h"
#include "lambda.h"
#include "io.h"
#include "probability.h"
#include "root_equilibrium_distribution.h"

reconstruction_process::reconstruction_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree,
    int max_family_size,
    int max_root_family_size, std::vector<int> rootdist,
    gene_family *gf,
    probability_calculator *p_calc,
    root_equilibrium_distribution* p_prior) : _gene_family(gf), _p_calc(p_calc), _p_prior(p_prior),
    process(ost, lambda, lambda_multiplier, p_tree, max_family_size, max_root_family_size, rootdist)
{
}

class child_multiplier
{
    std::map<clade *, std::vector<double> >& _L;
    double value;
    int _j;
public:
    child_multiplier(std::map<clade *, std::vector<double> >& L) : _L(L)
    {
        value = 1.0;
        _j = 0;
    }

    void set_index(int j)
    {
        _j = j;
    }

    void operator()(clade *c)
    {
        value *= _L[c][_j];
    }

    double result() const {
        return value;
    }
};

void reconstruction_process::reconstruct_leaf_node(clade * c, lambda * _lambda)
{
    auto& C = all_node_Cs[c];
    auto& L = all_node_Ls[c];
    C.resize(_max_family_size + 1);
    L.resize(_max_family_size + 1);

    double branch_length = c->get_parent()->get_branch_length();

    single_lambda * sl = dynamic_cast<single_lambda *>(_lambda);

    L.resize(_max_family_size + 1);

    int observed_count = _gene_family->get_species_size(c->get_taxon_name());
    fill(C.begin(), C.end(), observed_count);

    // i will be the parent size
    for (size_t i = 1; i < L.size(); ++i)
    {
        L[i] = _p_calc->get_from_parent_fam_size_to_c(sl->get_single_lambda(), branch_length, i, observed_count, NULL);
    }
}

void reconstruction_process::reconstruct_root_node(clade * c)
{
    auto& L = all_node_Ls[c];
    int C = -1;

    L.resize(_max_root_family_size + 1);

    child_multiplier cr(all_node_Ls);

    // i is the parent, j is the child
    for (size_t i = 1; i < L.size(); ++i)
    {
        double max_val = -1;

        for (size_t j = 1; j < L.size(); ++j)
        {
            cr.set_index(j);
            c->apply_to_descendants(cr);
            double val = cr.result() * _p_prior->compute(j);
            if (val > max_val)
            {
                max_val = val;
                C = j;
            }
        }

        L[i] = max_val;
    }

    cout << "C for tree is " << C << std::endl;
}

void reconstruction_process::reconstruct_internal_node(clade * c, lambda * _lambda)
{
    auto& C = all_node_Cs[c];
    auto& L = all_node_Ls[c];
    C.resize(_max_family_size + 1);
    L.resize(_max_family_size + 1);

    double branch_length = c->get_parent()->get_branch_length();

    single_lambda * sl = dynamic_cast<single_lambda *>(_lambda);

    L.resize(_max_family_size + 1);

    child_multiplier cr(all_node_Ls);

    // i is the parent, j is the child
    for (size_t i = 0; i < L.size(); ++i)
    {
        size_t max_j;
        double max_val = -1;
        if (c->is_root())
        {
            reconstruct_root_node(c);
        }
        else
        {
            for (size_t j = 0; j < L.size(); ++j)
            {
                cr.set_index(j);
                c->apply_to_descendants(cr);
                double val = cr.result() *_p_calc->get_from_parent_fam_size_to_c(sl->get_single_lambda(), branch_length, i, j, NULL);
                if (val > max_val)
                {
                    max_j = j;
                    max_val = val;
                }
            }

        }

        L[i] = max_val;
        C[i] = max_j;
    }

}


void reconstruction_process::operator()(clade *c)
{
    // all_node_Cs and all_node_Ls are hashtables where the keys are nodes and values are vectors of doubles


    single_lambda * sl = dynamic_cast<single_lambda *>(_lambda);
    if (!sl)
    {
        throw std::runtime_error("Cannot reconstruct with multiple lambdas yet");
    }

    if (c->is_leaf())
    {
        reconstruct_leaf_node(c, _lambda);
    }
    else if (c->is_root())
    {
        reconstruct_root_node(c);
    }
    else
    {
        reconstruct_internal_node(c, _lambda);
    }
}

void reconstruction_process::reconstruct()
{
    // Pupko's joint reconstruction algorithm
    _p_tree->apply_reverse_level_order(*this);
}

