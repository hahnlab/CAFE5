#ifndef GAMMA_H
#define GAMMA_H

#include <vector>

#define point_gamma(prob, alpha, beta) point_chi2(prob, 2.0*(alpha))/(2.0*(beta)) // Simon: what is this macro doing?

double point_chi2(double prob, double v);

double incomplete_gamma(double x, double alpha, double ln_gamma_alpha);

double point_normal(double prob);

int discrete_gamma(double freqK[], double rK[], double alpha, double beta, int K, int median);

void get_gamma(std::vector<double> &v_freq, std::vector<double> &v_rate, double alpha); // wrapper

#endif /* GAMMA_H */

