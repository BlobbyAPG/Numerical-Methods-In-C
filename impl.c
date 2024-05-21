#include<stdio.h>
#include"num_methods.h"
#include<math.h>


// Root-finding Methods Section:
// - Bisection Method
// - Newton-Raphson Method
// - Secant Method
// - Regula Falsi Method (False Position)
// - Fixed-Point Iteration Method
// - Halley's Method

double bisection(double (*f)(double), double a, double b, double tol, int max_iter){
    double fa, fb, p, fp;
    int i = 0;
    while(i < max_iter){
        p = (a + b) / 2;
        fa = f(a);
        fb = f(b);
        fp = f(p);
        if(fp == 0 || (b - a) / 2 < tol){
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

double newt_raph(double (*f)(double), double (*df)(double), double x0, double tol, int max_iter){
    double x = x0;
    int i = 0;
    while(i < max_iter){
        double fx = f(x);
        double dfx = df(x);
        if(fabs(dfx) < 1e-14){
            return x;
        }
        double dx = fx / dfx;
        x = x - dx;
        if(fabs(dx) < tol){
            return x;
        }
        i++;
    }
    return x;
}

double secant(double (*f)(double), double x0, double x1, double tol, int max_iter){
    double x = x1;
    double x_prev = x0;
    int i = 0;
    while(i < max_iter){
        double fx = f(x);
        double fx_prev = f(x_prev);
        if(fx - fx_prev == 0){
            return x;
        }
        double dx = fx * (x - x_prev) / (fx - fx_prev);
        x_prev = x;
        x = x - dx;
        if(fabs(dx) < tol){
            return x;
        }
        i++;
    }
    return x;
}

double reg_falsi(double (*f)(double), double a, double b, double tol, int max_iter){
    double fa, fb, p, fp;
    int i = 0;
    while(i < max_iter){
        fa = f(a);
        fb = f(b);
        p = (a * fb - b * fa) / (fb - fa);
        fp = f(p);
        if(fp == 0 || (b - a) / 2 < tol){
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

double fixed_pt_iter(double (*g)(double), double x0, double tol, int max_iter){
    double x = x0;
    int i = 0;
    while(i < max_iter){
        double gx = g(x);
        double dx = gx - x;
        x = gx;
        if(fabs(dx) < tol){
            return x;
        }
        i++;
    }
    return x;
}

double halleys(double (*f)(double), double (*df)(double), double (*ddf)(double), double x0, double tol, int max_iter){
    double x = x0;
    int i = 0;
    while(i < max_iter){
        double fx = f(x);
        double dfx = df(x);
        double ddfx = ddf(x);
        if(fabs(dfx) < 1e-14){
            return x;
        }
        double dx = 2 * fx * dfx / (2 * pow(dfx, 2) - fx * ddfx);
        x = x - dx;
        if(fabs(dx) < tol){
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
            if (j != i) {
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
            coefficients[i] = (coefficients[i] - coefficients[i - 1]) / (x[i] - x[i - j]);
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
        term *= (u - i + 1) / i;
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
        term *= (u + i - 1) / i;
        result += term * y[n - 1 - i];
    }
    return result;
}

// Linear Algebra Methods Section:
// - Jacobi Iterative Method
// - Gauss-Seidel Iterative Method

void jacobi(double **A, double *b, double *x, int n, double tol, int max_iter){
    double *x_new = malloc(n * sizeof(double));
    int i = 0;
    while(i < max_iter){
        for (int j = 0; j < n; j++) {
            x_new[j] = b[j];
            for (int k = 0; k < n; k++) {
                if (j != k) {
                    x_new[j] -= A[j][k] * x[k];
                }
            }
            x_new[j] /= A[j][j];
        }
        double error = 0;
        for (int j = 0; j < n; j++) {
            error += fabs(x_new[j] - x[j]);
            x[j] = x_new[j];
        }
        if (error < tol) {
            break;
        }
        i++;
    }
    free(x_new);
}

void gau_sei(double **A, double *b, double *x, int n, double tol, int max_iter){
    int i = 0;
    while(i < max_iter){
        for (int j = 0; j < n; j++) {
            double sum = b[j];
            for (int k = 0; k < n; k++) {
                if (j != k) {
                    sum -= A[j][k] * x[k];
                }
            }
            x[j] = sum / A[j][j];
        }
        double error = 0;
        for (int j = 0; j < n; j++) {
            error += fabs(A[j][j] * x[j] - b[j]);
        }
        if (error < tol) {
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
    return (f(x + h) - f(x)) / h;
}

double two_pt_bckd_diff(double (*f)(double), double x, double h){
    return (f(x) - f(x - h)) / h;
}

double three_pt_fwd_diff(double (*f)(double), double x, double h){
    return (-3 * f(x) + 4 * f(x + h) - f(x + 2 * h)) / (2 * h);
}

double three_pt_bckd_diff(double (*f)(double), double x, double h){
    return (3 * f(x) - 4 * f(x - h) + f(x - 2 * h)) / (2 * h);
}

double three_pt_cent_diff(double (*f)(double), double x, double h){
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
    return h * sum;
}

double simps(double (*f)(double), double a, double b, int n){
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        sum += 2 * f(a + i * h) * (i % 2 == 0 ? 2 : 4);
    }
    return h * sum / 3;
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






int main(){

}