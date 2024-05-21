#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

// Root-finding methods
double bisection(double (*f)(double), double a, double b, double tol, int max_iter);
double newton_raphson(double (*f)(double), double (*df)(double), double x0, double tol, int max_iter);
double secant(double (*f)(double), double x0, double x1, double tol, int max_iter);
double false_position(double (*f)(double), double a, double b, double tol, int max_iter);
double fixed_point_iteration(double (*g)(double), double x0, double tol, int max_iter);
double halleys_method(double (*f)(double), double (*df)(double), double (*ddf)(double), double x0, double tol, int max_iter);

// Interpolation and approximation
double lagrange_interpolation(double *x, double *y, int n, double x_interp);
double divided_differences(double *x, double *y, double *coefficients, int n);
double forward_interpolation(double *x, double *y, int n, double x_interp);
double backward_interpolation(double *x, double *y, int n, double x_interp);

// Linear algebra
void jacobi_method(double **A, double *b, double *x, int n, double tol, int max_iter);
void gauss_seidel_method(double **A, double *b, double *x, int n, double tol, int max_iter);

// Numerical differentiation
double two_point_forward_diff(double (*f)(double), double x, double h);
double two_point_backward_diff(double (*f)(double), double x, double h);
double three_point_forward_diff(double (*f)(double), double x, double h);
double three_point_backward_diff(double (*f)(double), double x, double h);
double three_point_central_diff(double (*f)(double), double x, double h);

// Numerical integration
double trapezoidal_rule(double (*f)(double), double a, double b, int n);
double simpsons_rule(double (*f)(double), double a, double b, int n);

// Ordinary Differential Equations (ODEs)
void euler_method(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void modified_euler_method(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void runge_kutta_2(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void runge_kutta_3(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);
void runge_kutta_4(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y);

// Aitken's Delta-Squared process
double aitkens_delta_squared(double *x, int n);

#endif // NUMERICAL_METHODS_H
