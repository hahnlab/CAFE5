#ifndef root_equilibrium_distribution_h
#define root_equilibrium_distribution_h

#include <vector>

class clade;
class gene_family;
class input_parameters;

class root_equilibrium_distribution
{
public:
    virtual float compute(int val) const = 0;
    virtual void initialize(std::vector<int> rootdist_vec) = 0;
};

class uniform_distribution : public root_equilibrium_distribution
{
    std::vector<int> _rootdist_vec; // in case the user wants to use a specific root size distribution for all simulations
public:
    virtual void initialize(std::vector<int> rootdist_vec)
    {
        _rootdist_vec = rootdist_vec;
    }

    virtual float compute(int val) const;   // creates uniform
};

class poisson_distribution : public root_equilibrium_distribution
{
    std::vector<double> poisson;
    double _poisson_lambda;
public:
    poisson_distribution(clade *p_tree, std::vector<gene_family> *p_gene_families);

    poisson_distribution(double poisson_lambda) : _poisson_lambda(poisson_lambda)
    {
    }
    virtual void initialize(std::vector<int> rootdist_vec);

    virtual float compute(int val) const
    {
        return poisson[val];
    }

};

root_equilibrium_distribution* root_eq_dist_factory(const input_parameters& my_input_parameters, clade *p_tree, std::vector<gene_family> *p_gene_families);
#endif
