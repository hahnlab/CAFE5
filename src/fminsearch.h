#ifndef FMINSEARCH_H
#define FMINSEARCH_H

typedef double(*math_func)(double* x, void* args);

struct FMinSearch
{
  int maxiters;
  int bymax;
  double rho, chi, psi, sigma;
  double tolx, tolf;
  double delta, zero_delta;

  int 	N, N1;
  int 	iters;
  double** v;
  double* fv;
  double** vsort;
  double* x_mean;
  double* x_r;
  double* x_tmp;
  int*    idx;

  void* args;
  math_func eq;
};

FMinSearch* fminsearch_new();
FMinSearch* fminsearch_new_with_eq(math_func eq, int Xsize, void* args);
void fminsearch_set_equation(FMinSearch* pfm, math_func eq, int Xsize, void* args);
int fminsearch_min(FMinSearch* pfm, double* X0);
double* fminsearch_get_minX(FMinSearch* pfm);
double fminsearch_get_minF(FMinSearch* pfm);

#endif
