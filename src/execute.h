#ifndef EXECUTE_H
#define EXECUTE_H

#include "io.h" 

struct input_parameters;

class execute {
public:
    //! Read in gene family data (-i)
    std::vector<gene_family> * read_gene_family_data(const input_parameters &my_parameters, int &max_family_size, int &max_root_family_size); 
};


#endif /* EXECUTE_H */
