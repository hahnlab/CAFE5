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
};


class matrix_cache_key {
    double _lambda;
    double _branch_length;
    size_t _size;
public:
    matrix_cache_key(int size, double some_lambda, double some_branch_length) :
        _size(size),
        _lambda(some_lambda),
        _branch_length(some_branch_length) { }
    bool operator<(const matrix_cache_key &o) const {
        return std::tie(_size, _branch_length, _lambda) < std::tie(o._size, o._branch_length, o._lambda);
    }
    double lambda() const {
        return _lambda;
    }
    double branch_length() const {
        return _branch_length;
    }
};


//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
Contains a map (_cache) that serves as a hash table to store precalculated values.
If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class matrix_cache {
private:
    std::map<matrix_cache_key, matrix*> _matrix_cache; //!< nested map that stores transition probabilities for a given lambda and branch_length (outer), then for a given parent and child size (inner)
public:
    double get_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int child_size) const;
    matrix get_matrix(int size, double branch_length, double lambda);
    void precalculate_matrices(int size, lambda* lambda, const std::set<double>& branch_lengths);

    int get_cache_size() const {
        return _matrix_cache.size();
    }

    ~matrix_cache();
};
#endif
