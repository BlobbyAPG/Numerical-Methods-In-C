#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

// Root-finding methods
double bisection(double (*f)(double), double a, double b, double tol, int max_iter);
double newt_raph(double (*f)(double), double (*df)(double), double x0, double tol, int max_iter);
double secant(double (*f)(double), double x0, double x1, double tol, int max_iter);
double reg_falsi(double (*f)(double), double a, double b, double tol, int max_iter);
double fixed_pt_iter(double (*g)(double), double x0, double tol, int max_iter);
double halleys(double (*f)(double), double (*df)(double), double (*ddf)(double), double x0, double tol, int max_iter);

// Interpolation and approximation
double lagrange(double *x, double *y, int n, double x_interp);
double divi_diffs(double *x, double *y, double *coefficients, int n);
double fwd_intpol(double *x, double *y, int n, double x_interp);
double bckwd_intpol(double *x, double *y, int n, double x_interp);

// Linear algebra
void jacobi(double **A, double *b, double *x, int n, double tol, int max_iter);
void gau_sei(double **A, double *b, double *x, int n, double tol, int max_iter);

// Numerical differentiation
double two_pt_forward_diff(double (*f)(double), double x, double h);
double two_pt_bckd_diff(double (*f)(double), double x, double h);
double three_pt_fwd_diff(double (*f)(double), double x, double h);
double three_pt_bckd_diff(double (*f)(double), double x, double h);
double three_pt_cent_diff(double (*f)(double), double x, double h);

// Numerical integration
double trap(double (*f)(double), double a, double b, int n);
double simps(double (*f)(double), double a, double b, int n);

// Ordinary Differential Equations (ODEs)
void euler(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void mod_euler(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void rk2(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void rk3(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void rk4(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);

// Aitken's Delta-Squared process
double atkn_dlta_sqd(double *x, int n);

#endif // NUMERICAL_METHODS_H
