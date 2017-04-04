#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

double unifrnd();
double the_probability_of_going_from_parent_fam_size_to_c(double lambda, int branch_length, int parent_size, int size);
double gammaln(double a);
double chooseln(double n, double k);

#endif
