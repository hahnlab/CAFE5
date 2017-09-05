#ifndef EXECUTE_H
#define EXECUTE_H

#include "io.h" 

struct input_parameters;

class execute {
public:
    //! Read in gene family data (-i)
    std::vector<gene_family> * read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size);
    
    //! Read in phylogenetic tree data (-t)
    clade * read_input_tree(const input_parameters &my_input_parameters, bool read_lambdas);
};


#endif /* EXECUTE_H */
