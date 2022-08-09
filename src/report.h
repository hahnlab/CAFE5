#ifndef REPORT_H
#define REPORT_H

#include <iosfwd>

class clade;
class lambda;
class gene_family;
class reconstruction;

template<typename T>
using clademap = std::map<const clade*, T>;

struct Report {
    const clade* p_tree = nullptr;
    const clade* p_lambda_tree = nullptr;
    const lambda* p_lambda = nullptr;
    clademap<float> average_expansion;

    struct delta_counts { int expansion; int decrease; int remain; };
    clademap<delta_counts> delta_count;

    void compute_expansion(const std::vector<gene_family>& gene_families, const reconstruction& reconstruct);
};

std::ostream& operator<<(std::ostream& ost, const Report& report);

#endif