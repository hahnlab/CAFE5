#include <algorithm>
#include <cmath>
#include <iostream>
#include "probability.h"

/* Useful links
1) http://www.rskey.org/gamma.htm # explanation for lgamma
2) http://www.physics.unlv.edu/~pang/cp_c.html # c code
*/

/* Necessary for old C implementation of what now is lgamma
#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
*/

/* Necessary for old C implementation of what now is lgamma 
double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
 24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
 -5.395239384953e-6 };
*/

/* Old C implementation of what now is lgamma */
/*
double gammaln(double a)
{
	int n;
	double p = __Qs[0];
	double a_add_5p5 = a + 5.5;
	for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
	return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}
*/

/* Old C implementation necessary for set_node_familysize_random. Now using uniform_real_distribution()
double unifrnd()
{
  return rand() / (RAND_MAX + 1.0); // rand() returns an int from 0 to RAND_MAX (which is defined in std); the +1.0 is there probably so that we do not draw exactly 1.
}
*/

double chooseln(double n, double r)
{

  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);
  return lgamma(n + 1) - lgamma(r + 1) - lgamma(n - r + 1);
}

/* Eqn. (1) in 2005 paper. Assumes u = lambda */
double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff)
{

  int m = std::min(c, s);
  double lastterm = 1;
  double p = 0.0;
  int s_add_c = s + c;
  int s_add_c_sub_1 = s_add_c - 1;
  int s_sub_1 = s - 1;

  for (int j = 0; j <= m; j++) {
    double t = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
    p += (exp(t) * lastterm); // Note that t is in log scale, therefore we need to do exp(t) to match Eqn. (1)
    lastterm *= coeff; // equivalent of ^j in Eqn. (1)
  }

  return std::max(std::min(p, 1.0), 0.0);
}

double the_probability_of_going_from_parent_fam_size_to_c(double lambda, double branch_length, int parent_size, int size)
{

  double alpha = lambda*branch_length / (1 + lambda*branch_length);
  double coeff = 1 - 2 * alpha;
  
  return birthdeath_rate_with_log_alpha(parent_size, size, log(alpha), coeff);
}
