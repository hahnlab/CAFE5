#ifndef MATRIX_CACHE_H
#define MATRIX_CACHE_H

#include <map>
#include <vector>
#include <set>

#include <assert.h>

class lambda;
class readwritelock;

class matrix
{
    std::vector<double> values;
    int _size;
public:
    matrix(int sz) : _size(sz)
    {
        values.resize(_size*_size);
    }
    void set(int x, int y, double val)
    {
        assert(x < _size);
        assert(y < _size);
        values[x*_size + y] = val;
    }
    double get(int x, int y) const
    {
        assert(x < _size);
        assert(y < _size);
        return values[x*_size + y];
    }
    int size() const {
        return _size;
    }
    bool is_zero() const;
    std::vector<double> multiply(const std::vector<double>& v, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size) const;
};


class matrix_cache_key {
    size_t _size;
    long _lambda;
    long _branch_length;
public:
    matrix_cache_key(int size, double some_lambda, double some_branch_length) :
        _size(size),
        _lambda(long(some_lambda * 1000000000)),    // keep 9 significant digits
        _branch_length(long(some_branch_length * 1000)) {} // keep 3 significant digits

    bool operator<(const matrix_cache_key &o) const {
        return std::tie(_size, _branch_length, _lambda) < std::tie(o._size, o._branch_length, o._lambda);
    }
    double lambda() const {
        return double(_lambda) / 1000000000.0;
    }
    double branch_length() const {
        return double(_branch_length) / 1000.0;
    }
};

std::vector<double> get_lambda_values(const lambda *p_lambda);

//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
Contains a map (_cache) that serves as a hash table to store precalculated values.
If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class matrix_cache {
private:
    std::map<matrix_cache_key, matrix*> _matrix_cache; //!< nested map that stores transition probabilities for a given lambda and branch_length (outer), then for a given parent and child size (inner)
    int _matrix_size;
public:
    double get_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int child_size) const;
    const matrix* get_matrix(double branch_length, double lambda) const;
    void precalculate_matrices(const std::vector<double>& lambdas, const std::set<double>& branch_lengths);

    int get_cache_size() const {
        return _matrix_cache.size();
    }

    int get_matrix_size() const {
        return _matrix_size;
    }

    void warn_on_saturation(std::ostream& ost);

    static bool is_saturated(double branch_length, double lambda);

    matrix_cache(int matrix_size) : _matrix_size(matrix_size) {}
    ~matrix_cache();
};
#endif
