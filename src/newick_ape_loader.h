#ifndef NEWICK_APE_LOADER_H
#define NEWICK_APE_LOADER_H

#include <vector>

class clade;
using cladevector = std::vector<const clade*>;

cladevector get_ape_order(const clade* p_tree);

#endif

