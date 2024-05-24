#include<stdio.h>
#include"num_methods.h"
#include<math.h>
#include<stdlib.h>


// Root-finding Methods Section:
// - Bisection Method
// - Newton-Raphson Method
// - Secant Method
// - Regula Falsi Method (False Position)
// - Fixed-Point Iteration Method
// - Halley's Method

double bisection(double (*f)(double), double a, double b, double err, int iter){
    double fa, fb, p, fp;
    int i = 0;
    while(i < iter){
        p = (a + b) / 2;
        fa = f(a);
        fb = f(b);
        fp = f(p);
        if(fp == 0 || (b - a) / 2 < err){
            return p;
        }
        if(fa * fp > 0){
            a = p;
        }else{
            b = p;
        }
        i++;
    }
    return p;
}

double newt_raph(double (*f)(double), double (*df)(double), double x0, double err, int iter){
    double x = x0;
    int i = 0;
    while(i < iter){
        double fx = f(x);
        double dfx = df(x);
        if(fabs(dfx) < 1e-14){
            return x;
        }
        double dx = fx / dfx;
        x = x - dx;
        if(fabs(dx) < err){
            return x;
        }
        i++;
    }
    return x;
}

double secant(double (*f)(double), double x0, double x1, double err, int iter){
    double x = x1;
    double x_prev = x0;
    int i = 0;
    while(i < iter){
        double fx = f(x);
        double fx_prev = f(x_prev);
        if(fx - fx_prev == 0){
            return x;
        }
        double dx = fx * (x - x_prev) / (fx - fx_prev);
        x_prev = x;
        x = x - dx;
        if(fabs(dx) < err){
            return x;
        }
        i++;
    }
    return x;
}

double reg_falsi(double (*f)(double), double a, double b, double err, int iter){
    double fa, fb, p, fp;
    int i = 0;
    while(i < iter){
        fa = f(a);
        fb = f(b);
        p = (a * fb - b * fa) / (fb - fa);
        fp = f(p);
        if(fp == 0 || (b - a) / 2 < err){
            return p;
        }
        if(fa * fp > 0){
            a = p;
        }else{
            b = p;
        }
        i++;
    }
    return p;
}

double fixed_pt_iter(double (*g)(double), double x0, double err, int iter){
    double x = x0;
    int i = 0;
    while(i < iter){
        double gx = g(x);
        double dx = gx - x;
        x = gx;
        if(fabs(dx) < err){
            return x;
        }
        i++;
    }
    return x;
}

double halleys(double (*f)(double), double (*df)(double), double (*ddf)(double), double x0, double err, int iter){
    double x = x0;
    int i = 0;
    while(i < iter){
        double fx = f(x);
        double dfx = df(x);
        double ddfx = ddf(x);
        if(fabs(dfx) < 1e-14){
            return x;
        }
        double dx = 2 * fx * dfx / (2 * pow(dfx, 2) - fx * ddfx);
        x = x - dx;
        if(fabs(dx) < err){
            return x;
        }
        i++;
    }
    return x;
}

// Interpolation & Approximation Methods Section:
// - Lagrange Interpolation
// - Divided Differences Interpolation
// - Forward Differences Interpolation
// - Backward Differences Interpoaltion

double lagrange(double *x, double *y, int n, double x_interp){
    double result = 0;
    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (j != i && x[i] != x[j]) { // Avoid division by zero
                term *= (x_interp - x[j]) / (x[i] - x[j]);
            }
        }
        printf("Term %d: %.4f\n", i+1, term); // Print intermediate result
        result += term;
    }
    return result;
}

double divi_diffs(double *x, double *y, double *coefficients, int n){
    for (int i = 0; i < n; i++) {
        coefficients[i] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            if (x[i] != x[i - j]) { // Avoid division by zero
                coefficients[i] = (coefficients[i] - coefficients[i - 1]) / (x[i] - x[i - j]);
            } else {
                // Handle repeated x values
                coefficients[i] = 0; // Set coefficient to zero
            }
        }
    }
    return coefficients[n - 1];
}

double fwd_intpol(double *x, double *y, int n, double x_interp){
    double h = x[1] - x[0];
    double result = y[0];
    double term = 1;
    double u = (x_interp - x[0]) / h;
    for (int i = 1; i < n; i++) {
        if (h != 0) { // Avoid division by zero
            term *= (u - i + 1) / i;
        } else {
            term = 0; // Set term to zero
        }
        result += term * y[i];
    }
    return result;
}

double bckwd_intpol(double *x, double *y, int n, double x_interp){
    double h = x[1] - x[0];
    double result = y[n - 1];
    double term = 1;
    double u = (x_interp - x[n - 1]) / h;
    for (int i = 1; i < n; i++) {
        if (h != 0) { // Avoid division by zero
            term *= (u + i - 1) / i;
        } else {
            term = 0; // Set term to zero
        }
        result += term * y[n - 1 - i];
    }
    return result;
}

// Linear Algebra Methods Section:
// - Jacobi Iterative Method
// - Gauss-Seidel Iterative Method

void jacobi(double **A, double *b, double *x, int n, double err, int iter){
    if (A == NULL || b == NULL || x == NULL) {
        // Handle null pointers
        return;
    }
    double *x_new = malloc(n * sizeof(double));
    if (x_new == NULL) {
        // Handle memory allocation failure
        return;
    }
    int i = 0;
    while(i < iter){
        for (int j = 0; j < n; j++) {
            x_new[j] = b[j];
            for (int k = 0; k < n; k++) {
                if (j != k) {
                    x_new[j] -= A[j][k] * x[k];
                }
            }
            if (A[j][j] == 0) {
                // Handle division by zero
                free(x_new);
                return;
            }
            x_new[j] /= A[j][j];
        }
        double error = 0;
        for (int j = 0; j < n; j++) {
            error += fabs(x_new[j] - x[j]);
            x[j] = x_new[j];
        }
        if (error < err) {
            break;
        }
        i++;
    }
    free(x_new);
}

void gau_sei(double **A, double *b, double *x, int n, double err, int iter){
    if (A == NULL || b == NULL || x == NULL) {
        // Handle null pointers
        return;
    }
    int i = 0;
    while(i < iter){
        for (int j = 0; j < n; j++) {
            double sum = b[j];
            for (int k = 0; k < n; k++) {
                if (j != k) {
                    sum -= A[j][k] * x[k];
                }
            }
            if (A[j][j] == 0) {
                // Handle division by zero
                return;
            }
            x[j] = sum / A[j][j];
        }
        double error = 0;
        for (int j = 0; j < n; j++) {
            error += fabs(A[j][j] * x[j] - b[j]);
        }
        if (error < err) {
            break;
        }
        i++;
    }
}

// Numerical Differentiation Methods Section:
// - Two Point Forward Difference Method
// - Two Point Backward Difference Method
// - Three Point Forward Difference Method
// - Three Point Backward Difference Method
// - Three Point Central Difference Method

double two_pt_fwd_diff(double (*f)(double), double x, double h){
    // Check if h is zero to avoid division by zero
    if (fabs(h) < 1e-10) {
        // Handle division by zero
        return NAN; // Return NaN (Not-a-Number) to indicate error
    }
    // Compute the two-point forward difference formula
    return (f(x + h) - f(x)) / h;
}

double two_pt_bckd_diff(double (*f)(double), double x, double h){
    // Check if h is zero to avoid division by zero
    if (fabs(h) < 1e-10) {
        // Handle division by zero
        return NAN; // Return NaN (Not-a-Number) to indicate error
    }
    // Compute the two-point backward difference formula
    return (f(x) - f(x - h)) / h;
}

double three_pt_fwd_diff(double (*f)(double), double x, double h){
    // Check if h is zero to avoid division by zero
    if (fabs(h) < 1e-10) {
        // Handle division by zero
        return NAN; // Return NaN (Not-a-Number) to indicate error
    }
    // Compute the three-point forward difference formula
    return (-3 * f(x) + 4 * f(x + h) - f(x + 2 * h)) / (2 * h);
}

double three_pt_bckd_diff(double (*f)(double), double x, double h){
    // Check if h is zero to avoid division by zero
    if (fabs(h) < 1e-10) {
        // Handle division by zero
        return NAN; // Return NaN (Not-a-Number) to indicate error
    }
    // Compute the three-point backward difference formula
    return (3 * f(x) - 4 * f(x - h) + f(x - 2 * h)) / (2 * h);
}

double three_pt_cent_diff(double (*f)(double), double x, double h){
    // Check if h is zero to avoid division by zero
    if (fabs(h) < 1e-10) {
        // Handle division by zero
        return NAN; // Return NaN (Not-a-Number) to indicate error
    }
    // Compute the three-point central difference formula
    return (f(x + h) - f(x - h)) / (2 * h);
}


// Numerical Integration Methods Section:
// - Trapezoidal Rule
// - Simpson's Rule

double trap(double (*f)(double), double a, double b, int n){
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++) {
        sum += f(a + i * h);
    }
    int result = 0;
    return result = h * sum;

    // Calculate the true value of the integral
    double true_value = 1.0 / 3.0 * (pow(b, 3) - pow(a, 3));

    // Calculate the absolute error
    double abs_error = fabs(true_value - result);

    // Calculate the relative error
    double rel_error = fabs(abs_error / true_value);

    printf("Absolute Error: %lf\n", abs_error);
    printf("Relative Error: %lf\n", rel_error);

    return result;
}

double simps(double (*f)(double), double a, double b, int n){
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        sum += 2 * f(a + i * h) * (i % 2 == 0 ? 2 : 4);
    }
    int result = 0;
    return result = h * sum / 3;

// Calculate the true value of the integral
    double true_value = 1.0 / 3.0 * (pow(b, 3) - pow(a, 3));

    // Calculate the absolute error
    double abs_error = fabs(true_value - result);

    // Calculate the relative error
    double rel_error = fabs(abs_error / true_value);

    printf("Absolute Error: %lf\n", abs_error);
    printf("Relative Error: %lf\n", rel_error);
}

// Ordinary Differential Equations Section:
// - Euler's Method
// - Modified Euler's Method
// - Runge-Kutta's 2nd Order
// - Runge-Kutta's 3rd Order
// - Runge-Kutta's 4th Order

void euler(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y){
    x[0] = x0;
    y[0] = y0;
    for (int i = 1; i <= n; i++) {
        x[i] = x[i - 1] + h;
        y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1]);
    }
}

void mod_euler(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y){
    x[0] = x0;
    y[0] = y0;
    for (int i = 1; i <= n; i++) {
        x[i] = x[i - 1] + h;
        double k1 = h * f(x[i - 1], y[i - 1]);
        double k2 = h * f(x[i], y[i - 1] + k1);
        y[i] = y[i - 1] + 0.5 * (k1 + k2);
    }
}

void rk2(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y){
    x[0] = x0;
    y[0] = y0;
    for (int i = 1; i <= n; i++) {
        x[i] = x[i - 1] + h;
        double k1 = h * f(x[i - 1], y[i - 1]);
        double k2 = h * f(x[i], y[i - 1] + k1);
        y[i] = y[i - 1] + 0.5 * (k1 + k2);
    }
}

void rk3(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y){
    x[0] = x0;
    y[0] = y0;
    for (int i = 1; i <= n; i++) {
        x[i] = x[i - 1] + h;
        double k1 = h * f(x[i - 1], y[i - 1]);
        double k2 = h * f(x[i - 1] + h / 2, y[i - 1] + k1 / 2);
        double k3 = h * f(x[i], y[i - 1] - k1 + 2 * k2);
        y[i] = y[i - 1] + (k1 + 4 * k2 + k3) / 6;
    }
}

void rk4(double (*f)(double, double), double x0, double y0, double h, int n, double *x, double *y){
    x[0] = x0;
    y[0] = y0;
    for (int i = 1; i <= n; i++) {
        x[i] = x[i - 1] + h;
        double k1 = h * f(x[i - 1], y[i - 1]);
        double k2 = h * f(x[i - 1] + h / 2, y[i - 1] + k1 / 2);
        double k3 = h * f(x[i - 1] + h / 2, y[i - 1] + k2 / 2);
        double k4 = h * f(x[i], y[i - 1] + k3);
        y[i] = y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
}

// Aitken's Delta-Squared Acceleration Process
double atkn_dlta_sqd(double *x, int n){
    double *y = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        y[i] = x[i];
    }
    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < n - i - 2; j++) {
            y[j] = y[j] - pow(y[j + 1] - y[j], 2) / (y[j + 2] - 2 * y[j + 1] + y[j]);
        }
    }
    double result = y[0];
    free(y);
    return result;
}


int main(){
    // Root-Finding Methods
    printf("Root-Finding Methods:\n");

    // Bisection Method
    printf("Bisection Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_bisection = bisection(f, 1, 2, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_bisection);

    // Newton-Raphson Method
    printf("Newton-Raphson Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_newton_raphson = newt_raph(f, df, 1.5, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_newton_raphson);

    // Secant Method
    printf("Secant Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_secant = secant(f, 1.5, 1.6, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_secant);

    // Regula Falsi Method
    printf("Regula Falsi Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_reg_falsi = reg_falsi(f, 1, 2, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_reg_falsi);

    // Fixed-Point Iteration Method
    printf("Fixed-Point Iteration Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_fixed_pt_iter = fixed_pt_iter(g, 1.5, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_fixed_pt_iter);

    // Halley's Method
    printf("Halley's Method:\n");
    printf("Iteration\tRoot\tConvergence Rate\n");
    double root_halleys = halleys(f, df, ddf, 1.5, 1e-6, 100);
    printf("%d\t\t%.6lf\n", iter, root_halleys);

    // Numerical Differentiation Methods
    printf("Numerical Differentiation Methods:\n");

    // Two-Point Forward Difference Method
    printf("Two-Point Forward Difference Method:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double result_two_pt_fwd_diff = two_pt_fwd_diff(f, M_PI / 4, 0.1);
    printf("%d\t\t%.6lf\n", iter, result_two_pt_fwd_diff);

    // Two-Point Backward Difference Method
    printf("Two-Point Backward Difference Method:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double result_two_pt_bckd_diff = two_pt_bckd_diff(f, M_PI / 4, 0.1);
    printf("%d\t\t%.6lf\n", iter, result_two_pt_bckd_diff);

    // Numerical Integration Methods
    printf("Numerical Integration Methods:\n");

    // Trapezoidal Rule
    printf("Trapezoidal Rule:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double result_trap = trap(f_sq, 0, 1, 10);
    printf("%d\t\t%.6lf\n", iter, result_trap);

    // Simpson's Rule
    printf("Simpson's Rule:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double result_simps = simps(f_sq, 0, 1, 10);
    printf("%d\t\t%.6lf\n", iter, result_simps);

    // Ordinary Differential Equations (ODEs)
    printf("Ordinary Differential Equations (ODEs):\n");

    double result;
// Example functions for testing
    double f(double x) {
        return x * x - 4; // x^2 - 4
    }

    double df(double x) {
        return 2 * x; // Derivative of x^2 - 4: 2x
    }

    double ddf(double x) {
        return 2; // Derivative of x^2 - 4: 2x
    }

    double g(double x) {
        return x * x * x - 2 * x - 5; // x^3 - 2x - 5
    }

    double g_diff(double x) {
        return 3 * x * x - 2; // Derivative of x^3 - 2x - 5: 3x^2 - 2
    }

// Example ordinary differential equation
    double f_diff(double x, double y) {
        return x + y; // Differential equation: y' = x + y
    }

    double f_sq(double x) {
        return x * x; // f(x) = x^2
    }

    int iter = 12;

    // Euler's Method
    printf("Euler's Method:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double x_euler[11], y_euler[11];
    euler(f_diff, 0, 1, 0.1, 10, x_euler, y_euler);
    for (int i = 0; i <= 10; i++) {
        printf("%d\t\t%.6lf\n", i, y_euler[i]);
    }

    // Modified Euler's Method
    printf("Modified Euler's Method:\n");
    printf("Iteration\tResult\tConvergence Rate\n");
    double x_mod_euler[11], y_mod_euler[11];
    mod_euler(f_diff, 0, 1, 0.1, 10, x_mod_euler, y_mod_euler);
    for (int i = 0; i <= 10; i++) {
        printf("%d\t\t%.6lf\n", i, y_mod_euler[i]);
    }
}