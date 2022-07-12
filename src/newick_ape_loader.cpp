/* Newick tree loader that reproduces the structure that */
/* the R Ape package uses */
/* based on work by Emmanuel Paradis and Klaus Schliep */

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

#include "doctest.h"

#include "clade.h"

using namespace std;


struct tree_data
{
	vector<int> edge;	// edge holds pairs of nodes that have edges: edge[0] and edge[edge.size()]/2 have an edge
	vector<double> edge_length;
	vector<string> node_label;
	vector<string> tip_label;
	vector<double> root_edge;
};

void extract_portion_Newick(string x, int a, int b, char* y)
{
	int i, j;

	for (i = a, j = 0; i <= b; i++, j++) y[j] = x[i];

	y[j] = '\0';
}

void decode_internal_edge(string x, int a, int b, char* lab, double* w)
{
	int co = a;
	char* endstr, str[100];

	while (x[co] != ':' && co <= b) co++;

	if (a == co) lab[0] = '\0'; /* if no node label */
	else extract_portion_Newick(x, a, co - 1, lab);
	if (co < b) {
		extract_portion_Newick(x, co + 1, b, str);
		*w = strtod(str, &endstr);
	}
	else *w = NAN;
}

void decode_terminal_edge(string x, int a, int b, char* tip, double* w)
{
	int co = a;
	char* endstr, str[100];

	while (x[co] != ':' && co <= b) co++;

	extract_portion_Newick(x, a, co - 1, tip);
	if (co < b) {
		extract_portion_Newick(x, co + 1, b, str);
		*w = strtod(str, &endstr);
	}
	else *w = NAN;
}

#define ADD_INTERNAL_EDGE            \
    phy.edge[j] = curnode;                  \
    phy.edge[j + nedge] = curnode = ++node; \
    stack_internal[k++] = j;         \
    j++

#define GO_DOWN                                                  \
    decode_internal_edge(x, ps + 1, pt - 1, lab, &tmpd);         \
    phy.node_label[curnode - 1 - ntip] = lab;		\
    l = stack_internal[--k];					 \
    phy.edge_length[l] = tmpd;                                                \
    curnode = phy.edge[l]

void add_terminal_edge_tiplabel(tree_data& phy, int& curtip, int& j, int curnode, string x, int pr, int ps, char* tip, int nedge)
{
	double tmpd;
	phy.edge[j] = curnode;
	decode_terminal_edge(x, pr + 1, ps - 1, tip, &tmpd);
	if (phy.tip_label.size() < (size_t)curtip) phy.tip_label.resize(curtip);
	phy.tip_label[curtip - 1] = tip;
	phy.edge[j + nedge] = curtip;
	phy.edge_length[j] = tmpd;
	curtip++;
	j++;
}

#define ADD_TERMINAL_EDGE_TIPLABEL add_terminal_edge_tiplabel(phy, curtip, j, curnode, x, pr, ps, tip, nedge)

tree_data treeBuild(string nwk)
{
	string x;
	int n, i, ntip = 1, nleft = 0, nright = 0, nedge, curnode, node, j, nsk = 0, ps, pr, pt, l, k, stack_internal[10000], curtip = 1;
	int nnode = 0;
	double tmpd;
	char lab[512], tip[512];
	tree_data phy;
	phy.tip_label.resize(1);

	x = nwk;
	n = x.size();
	vector<int> skeleton(n);
	for (i = 0; i < n; i++) {

		if (x[i] == '(') {

			skeleton[nsk] = i;
			nsk++;
			nleft++;
			continue;
		}
		if (x[i] == ',') {

			skeleton[nsk] = i;
			nsk++;
			ntip++;
			continue;
		}
		if (x[i] == ')') {

			skeleton[nsk] = i;
			nsk++;
			nright++;
			nnode++;
		}
	}
	if (nleft != nright) perror("numbers of left and right parentheses in Newick string not equaln");
	nedge = ntip + nnode - 1;
	phy.edge.resize(nedge * 2);
	phy.edge_length.resize(nedge);
	phy.node_label.resize(nnode);

	curnode = node = ntip + 1;	// curnode starts at number of tips + 1
	k = j = 0;

	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				ADD_TERMINAL_EDGE_TIPLABEL;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1];
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE_TIPLABEL;
				GO_DOWN;
				continue;
			}
			if (x[pr] == '(') {
				ADD_TERMINAL_EDGE_TIPLABEL;
				GO_DOWN;
				continue;
			}
			if (x[pr] == ')') {
				GO_DOWN;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE_TIPLABEL;
	}

	if (ps < n - 2) {
		i = ps + 1;
		while (i < n - 2 && x[i] != ':') i++;
		if (i < n - 2) {
			decode_internal_edge(x, ps + 1, n - 2, lab, &tmpd);
			phy.root_edge.push_back(tmpd);
			phy.node_label[0] = lab;
		}
		else {
			extract_portion_Newick(x, ps + 1, n - 2, lab);
			phy.node_label[0] = lab;
		}
	}

	return phy;
}

set<int> get_connections(tree_data& td, int node_id)
{
	set<int> result;
	size_t split = td.edge.size() / 2;
	for (size_t i = 0; i < td.edge.size(); ++i)
	{
		if (td.edge[i] == node_id)
		{
			if (i < split)
				result.insert(td.edge[i + split]);
			else
				result.insert(td.edge[i - split]);
		}
	}
	return result;
}

cladevector get_ape_order(const clade *p_tree)
{
	auto t = treeBuild(p_tree->get_source_newick());
	clademap<int> indices;
	indices[p_tree] = t.tip_label.size() + 1;
	for (size_t i = 0; i < t.tip_label.size(); ++i)
		indices[p_tree->find_descendant(t.tip_label[i])] = i+1;

	for (auto it = p_tree->reverse_level_begin(); it != p_tree->reverse_level_end(); ++it)
	{
		const clade* c = *it;
		if (c->is_leaf() || c->is_root()) continue;
		// find an edge of the child that isn't already in our list - that must be the child's parent, e.g., our node
		for (auto edge : get_connections(t, indices[*c->descendant_begin()]))
		{
			if (find_if(indices.begin(), indices.end(), [edge](clademap<int>::value_type p ) { return p.second == edge; }) == indices.end())
			{
				indices[c] = edge;
				break;
			}
		}
	}
	auto comp = [](std::pair<const clade* const, int>& a, std::pair<const clade* const, int>& b) { return a.second < b.second;  };
	int count = max_element(indices.begin(), indices.end(),comp)->second;
	cladevector result(count+1);
	for (auto& p : indices)
		result[p.second] = p.first;

	return result;
}

TEST_CASE("treeBuild")
{
	auto t = treeBuild(" ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)");

	CHECK_EQ(44, t.edge.size());
}

TEST_CASE("treeBuild orders tips as they appear in the newick string")
{
	auto t = treeBuild("((Parasteatoda_tepidariorum<285>_1:138.531,Stegodyphus_mimosarum<284>_1:138.531)<297>_1:187.093,Centruroides_sculpturatus<296>_3:325.625)<313>_1:46.2175;");

	REQUIRE_EQ(3, t.tip_label.size());
	CHECK_EQ("Parasteatoda_tepidariorum<285>_1", t.tip_label[0]);
	CHECK_EQ("Stegodyphus_mimosarum<284>_1", t.tip_label[1]);
	CHECK_EQ("Centruroides_sculpturatus<296>_3", t.tip_label[2]);
}

TEST_CASE("get_connections")
{
	auto t = treeBuild(" ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)");
	auto con = get_connections(t, 1);	// leaf
	CHECK_EQ(1, con.size());

	con = get_connections(t, 13);	    // root
	CHECK_EQ(2, con.size());

	con = get_connections(t, 14);	    // internal
	CHECK_EQ(3, con.size());
}

TEST_CASE("get_ape_order")
{
	unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
	auto actual = get_ape_order(p_tree.get());

	vector<string> nodes({ "A", "B", "C", "D", "ABCD", "AB", "CD" });
	cladevector expected(nodes.size()+1);
	transform(nodes.begin(), nodes.end(), expected.begin()+1, [&p_tree](string s) { return p_tree->find_descendant(s); });
	CHECK_EQ(expected, actual);

}