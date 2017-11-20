#ifndef CHISQUARE_H
#define CHISQUARE_H

#include <cmath>

#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
#define __EPS__		1e-8
/*
* URL : http://www.rskey.org/gamma.htm
*       http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/pdf.htm
*       http://www.physics.unlv.edu/~pang/cp_c.html
*/
extern double __Qs[];
/*
* URL: http://www.hpcsoft.com/products/MathSoL/specialFunction/gammaIn.html
* \int^{\infty}_{0} t^{a-1}e^{-t}\;dt
*/
inline double gammaln(double a)
{
    int n;
    double p = __Qs[0];
    double a_add_5p5 = a + 5.5;
    for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
    return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}

/*
* For x < \alpha + 1
*  gamma(x,\alpha) = e^{-x}x^{\alpha} \sum_{i=0}^{\infty} \frac{\Gamma(\alpha)}{\Gamma(\alpha+1+i)} x^{i}
*/

inline double incgammaln_lower(double x, double a)
{
    int i;
    double p = 1 / a;
    double t = 1 / a;
    for (i = 1; i < 1000; i++)
    {
        t *= x / (a + i);
        if (t < __EPS__) break;
        p += t;
    }
    return i == 1000 ? gammaln(a) : log(p) + a * log(x) - x;
}

/* \alpha > 0
* P(x,\alpha) = \frac{1}{\Gamma{\alpha}\int_0^{x} x^{\alpha-1}e^{-t}\;dt
*/
inline double gammaincln(double x, double a)
{
    return incgammaln_lower(x, a) - gammaln(a);
}

inline double gammainc(double x, double a)
{
    return exp(gammaincln(x, a));
}
/*
* \frac{ * \Gamma(\alpha,\frac{x}{\beta} } { \Gamma(\alpha) }
*/
inline double gamcdf(double x, double alpha, double beta)
{
    return gammainc(x / beta, alpha);
}

inline double chi2cdf(double x, int df)
{
    return gamcdf(x, df / 2.0, 2);
}
#endif