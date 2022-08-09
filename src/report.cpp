#include <ostream>
#include <iterator>
#include <numeric>

#define ELPP_NO_CHECK_MACROS
#include "easylogging++.h"
#include "doctest.h"

#include "clade.h"
#include "lambda.h"
#include "report.h"
#include "newick_ape_loader.h"
#include "gene_family.h"
#include "core.h"

using namespace std;

template<typename InputIt>
std::string join(InputIt first, InputIt last,
    const std::string& separator = ", ",
    const std::string& concluder = "")
{
    const std::string empty{};
    auto sep = &empty;

    std::stringstream ss;

    while (first != last)
    {
        ss << *sep << *first++;
        sep = &separator;
    }

    ss << concluder;

    return ss.str();
}

template<typename T, typename U>
std::ostream& operator<<(ostream& ost, const std::pair<T, U>& pair)
{
    ost << "(" << pair.first << "," << pair.second << ")";
    return ost;
}

std::ostream& operator<<(ostream& ost, const Report& report)
{
    //bool has_pvalues = false;
    //bool has_likelihoods = false;

    ost << "Tree:";
    report.p_tree->write_newick(ost, [](const clade* c)
        {
            ostringstream ost;
            ost << (c->is_leaf() ? c->get_taxon_name() : "") << ":" << c->get_branch_length();
            return ost.str();
        });
    ost << "\n";

    if (report.p_lambda)
    {
        ost << "Lambda:\t";
        auto vals = get_lambda_values(report.p_lambda);
        copy(vals.begin(), vals.end(), ostream_iterator<double>(ost, "\t"));
    }
    ost << "\n";
    if (report.p_lambda_tree)
    {
        ost << "Lambda tree:\t";
        report.p_lambda_tree->write_newick(ost, [](const clade* c)
            {
                ostringstream ost;
                ost << (c->is_leaf() ? c->get_taxon_name() : "") << ":" << c->get_branch_length();
                return ost.str();
            });
        ost << "\n";
    }

    auto t = get_ape_order(report.p_tree);
    map<const clade*, int> order;
    for (size_t i = 0; i < t.size(); ++i)
    {
        if (t[i]) order[t[i]] = i;
    }
    ost << "# IDs of nodes:";
    report.p_tree->write_newick(ost, [&order](const clade* c)
        {
            ostringstream ost;
            ost << (c->is_leaf() ? c->get_taxon_name() : "") << "<" << order[c] << ">";
            return ost.str();
        });
    ost << "\n";


    ost << "# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ";
    report.p_tree->apply_prefix_order([&ost, &order](const clade* c) {
        if (!c->is_leaf())
        {
            vector<int> v;
            transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [&order](const clade* c) { return order[c];  });
            ost << "(" << join(v.begin(), v.end(), ",") << ") ";
        }
    });

    ost << "\n";

    if (!report.average_expansion.empty())
    {
        ost << "Average Expansion:";
        report.p_tree->apply_prefix_order([&ost, &report](const clade* c) {
            if (!c->is_leaf())
            {
                vector<float> v;
                transform(c->descendant_begin(), c->descendant_end(), back_inserter(v), [&report](const clade* c) { return report.average_expansion.at(c);  });
                ost << "(" << join(v.begin(), v.end(), ",") << ") ";
            }
            });
        ost << "\n";
    }
#if 0

    write_viterbi(ost, report);

    write_families_header(ost, has_pvalues, has_likelihoods);

    copy(report.family_line_items.begin(), report.family_line_items.end(), ostream_iterator<family_line_item>(ost, "\n"));
#endif
    return ost;
}

void Report::compute_expansion(const std::vector<gene_family>& gene_families, const reconstruction& reconstruct)
{
    p_tree->apply_prefix_order([gene_families, &reconstruct, this](const clade* c) {
        if (!c->is_root())
        {
            int total = accumulate(gene_families.begin(), gene_families.end(), 0, [&reconstruct, c](int cur, const gene_family& gf) {
                return cur + reconstruct.get_difference_from_parent(gf, c);
                });
            this->average_expansion[c] = float(total) / float(gene_families.size());
        }
        });
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("Report writes tree")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Tree:((A:1,B:1):1,(C:1,D:1):1):0");
}

TEST_CASE("Report writes lambda tree")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
    unique_ptr<clade> p_tree2(parse_newick("((A:2,B:2):2,(C:3,D:3):1);"));

    Report r;
    r.p_tree = p_tree.get();
    r.p_lambda_tree = p_tree2.get();
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda tree:\t((A:2,B:2):2,(C:3,D:3):1):0");
}

TEST_CASE("Report writes lambdas")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    r.p_lambda = new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 });
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "Lambda:\t0.5\t0.3\t0.9");
}

TEST_CASE("Report writes IDs of nodes")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    r.p_lambda = new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 });
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, "# IDs of nodes:((A<1>,B<2>)<6>,(C<3>,D<4>)<7>)<5>");
}

TEST_CASE("Report writes node pairs")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    Report r;
    r.p_tree = p_tree.get();
    r.p_lambda = new multiple_lambda(map<std::string, int>(), { 0.5, 0.3, 0.9 });
    ostringstream ost;
    ost << r;

    CHECK_STREAM_CONTAINS(ost, " (node ID, node ID): (6,7) (1,2) (3,4)");
}

class mock_reconstruction : public reconstruction   
{
    virtual int get_node_count(const gene_family& gf, const clade* c) const override 
    { 
        return values.at(gf.id())->find_descendant(c->get_taxon_name())->get_branch_length();
    }
    map<string, const clade*> values;
public:
    void add(string id, string nwk) 
    {
        values[id] = parse_newick(nwk);
    }
};

TEST_CASE("Report compute_expansion")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    mock_reconstruction r2;
    r2.add("f1", "((A:1,B:3):2,(C:2,D:1):2);");
    r2.add("f2", "((A:1,B:5):2,(C:7,D:1):4);");
    r2.add("f3", "((A:1,B:5):2,(C:3,D:1):3);");
    r2.add("f4", "((A:1,B:5):2,(C:5,D:1):3);");
    Report r;
    r.p_tree = p_tree.get();

    vector<gene_family> f(4);
    f[0].set_id("f1");
    f[1].set_id("f2");
    f[2].set_id("f3");
    f[3].set_id("f4");

    r.compute_expansion(f, r2);
    CHECK_EQ(6, r.average_expansion.size());
    CHECK_EQ(-1, r.average_expansion[p_tree->find_descendant("A")]);
    CHECK_EQ(2.5, r.average_expansion[p_tree->find_descendant("B")]);
    CHECK_EQ(1.25, r.average_expansion[p_tree->find_descendant("C")]);
    CHECK_EQ(-2, r.average_expansion[p_tree->find_descendant("D")]);

    ostringstream ost;
    ost << r;
    CHECK_STREAM_CONTAINS(ost, "Average Expansion:(2,3) (-1,2.5) (1.25,-2)");
}

