#include <algorithm>
#include <cmath>

#include "probability.h"

#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)

/*
* URL : http://www.rskey.org/gamma.htm
*       http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/pdf.htm
*       http://www.physics.unlv.edu/~pang/cp_c.html
*/

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


double unifrnd()
{
	return rand() / (RAND_MAX + 1.0);
}

/*
* URL: http://www.hpcsoft.com/products/MathSoL/specialFunction/gammaIn.html
* \int^{\infty}_{0} t^{a-1}e^{-t}\;dt
*/
double gammaln(double a)
{
	int n;
	double p = __Qs[0];
	double a_add_5p5 = a + 5.5;
	for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
	return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}

double chooseln(double n, double r)
{
	if (r == 0 || (n == 0 && r == 0)) return 0;
	else if (n <= 0 || r <= 0) return log(0);
	return gammaln(n + 1) - gammaln(r + 1) - gammaln(n - r + 1);
}

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff)
{
	int m = std::min(c, s);

	double lastterm = 1;
	double p = 0.0;
	int s_add_c = s + c;
	int s_add_c_sub_1 = s_add_c - 1;
	int s_sub_1 = s - 1;
	for (int j = 0; j <= m; j++)
	{
		double t = chooseln(s, j) + chooseln(s_add_c_sub_1 - j, s_sub_1) + (s_add_c - 2 * j)*log_alpha;
		p += (exp(t) * lastterm);
		lastterm *= coeff;
	}

	return std::max(std::min(p, 1.0), 0.0);
}

double the_probability_of_going_from_parent_fam_size_to_c(double lambda, int branch_length, int parent_size, int size)
{
	double alpha = lambda*branch_length / (1 + lambda*branch_length);
	double coeff = 1 - 2 * alpha;

	birthdeath_rate_with_log_alpha(parent_size, size, log(alpha), coeff);
	return 0.05;
}
