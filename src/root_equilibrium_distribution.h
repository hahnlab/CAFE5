#ifndef root_equilibrium_distribution_h
#define root_equilibrium_distribution_h

#include <vector>
#include <map>

#include "gene_family.h"

class clade;
struct input_parameters;
class root_distribution;

class root_equilibrium_distribution
{
public:
    virtual float compute(size_t val) const = 0;
    virtual ~root_equilibrium_distribution() {}

    virtual int select_root_size(int family_number) const = 0;
};

class uniform_distribution : public root_equilibrium_distribution
{
    int _max_root_family_size = 0;
public:
    uniform_distribution(int max_root_family_size) : _max_root_family_size(max_root_family_size)
    {
    }

    virtual float compute(size_t val) const override
    {
        return 1.0 / _max_root_family_size;
    }

    int select_root_size(int family_number) const override;
};

class specified_distribution : public root_equilibrium_distribution
{
    std::vector<int> _vectorized_distribution; 
public:
    specified_distribution(const std::map<int, int>& root_distribution);

    virtual float compute(size_t val) const override;   // creates uniform

    int select_root_size(int family_number) const override;

    void resize(size_t new_size);
};

class poisson_distribution : public root_equilibrium_distribution
{
    std::vector<double> poisson;
    double _poisson_lambda;
public:
    poisson_distribution(std::vector<gene_family> *p_gene_families, int num_values);

    poisson_distribution(double poisson_lambda, int num_values);

    virtual float compute(size_t val) const override
    {
        if (val >= poisson.size())
            return 0;

        return poisson[val];
    }

    int select_root_size(int family_number) const override;
};

root_equilibrium_distribution* root_eq_dist_factory(const input_parameters& my_input_parameters, std::vector<gene_family> *p_gene_families, const std::map<int, int>& root_distribution, int max_root_family_size);
#endif
