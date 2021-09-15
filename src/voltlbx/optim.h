#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace voltlbx
{

    struct OptimBracket
    {
        double a;
        double b;
        double c;
        double fa;
        double fb;
        double fc;
    };


    template <typename F>
    double optimize_brent(F&& func, OptimBracket& br, double tol, double* fmin = 0)
    {
        static const double zeps = std::numeric_limits<double>::epsilon() / 1024;
        static const double cgold = 0.3819660;
        double a = std::min(br.a, br.c);
        double b = std::max(br.a, br.c);
        double v = br.b;
        double w = br.b;
        double x = br.b;
        double fv = br.fb;
        double fw = br.fb;
        double fx = br.fb;
        double d = 0;
        double e = 0;
        auto shift3 = [](double& a, double& b, double& c, double& d) {
            a = b;
            b = c;
            c = d;

        };
        for (unsigned j = 0; j < 100; j++) {
            double xm = 0.5 * (a + b);
            double toll = tol * std::abs(x) + zeps;
            double tol2 = 2 * toll;
            if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                if (fmin)
                    *fmin = fx;
                return x;
            }
            if (std::abs(e) > toll) {
                double r = (x - w) * (fx - fv);
                double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2 * (q - r);
                if (q > 0)
                    p = -p;
                q = std::abs(q);
                double etmp = e;
                e = d;
                if (std::abs(p) >= std::abs(0.5 * q * etmp) || p <= q * (a - x)
                    || p >= q * (b - q))
                {
                    d = cgold * (e = (x >= xm ? a - x : b - x));
                }
                else {
                    d = p / q;
                    double u = x + d;
                    if (u - a < tol2 || b - u < tol2)
                        d = std::copysign(toll, xm - x);
                }
            }
            else {
                d = cgold * (e = (x > xm ? a - x : b - x));
            }

            double u = std::abs(d) >= toll ? x + d : x + std::copysign(toll, d);
            double fu = func(u);
            if (fu <= fx) {
                if (u >= x)
                    a = x;
                else
                    b = x;
                shift3(v, w, x, u);
                shift3(fv, fw, fx, fu);

            }
            else {
                if (u < x)
                    a = u;
                else
                    b = u;
                if (fu <= fw || w == x) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }
        throw std::runtime_error("tooymany_ iterations");
    }


    class MinPackWrapper
    {
    public:

        template<typename F>
        MinPackWrapper(F&& f)
            : func(std::forward<F>(f))
        {}

        int solve(std::size_t n, double*x, std::size_t m, double *f, double tol);

    private:
        static int wrapper(void *, int m, int n, const double *x, double *f, int);

        std::exception_ptr error;
        std::function<void(std::size_t, const double *, std::size_t, double *)> func;
    };
    
    /* Optimize |f|^2 by Levenberg-Marquard algorithm.       
       https://github.com/devernay/cminpack
        Return code with followwing interpretation :
        INFO = 0 Improper input parameters.
        INFO = 1 Algorithm estimates that relative error in the sum of squares is at most TOL.
        INFO = 2 Algorithm estimates that relatives error between X and the solution is at most TOL.
        INFO = 3 Conditions for INFO = 1 nad INFO = 2 both hold.
        INFO = 4 FVEC is orthogonal to the Jacobian to machine precision.
        INFO = 5 Number of calls to FUNC has reached or exceeded 200*(N+1).
        INFO = 6 TOL is too small, no further reduction in sum of squares is possible.
        INFO = 7 TOL is too small, no further reduction in approximate solution X is possible.
    */
    template<typename F>
    int levmarq(F&& func, std::size_t n, double* x, std::size_t m, double *f, double tol)
    {
        auto w = MinPackWrapper(std::forward<F>(func));
        return w.solve(n, x, m, f, tol);
    }

    /* Optimize |f|^2 by Levenberg-Marquard algorithm.
       https://github.com/devernay/cminpack
        Return code with followwing interpretation :
        INFO = 0 Improper input parameters.
        INFO = 1 Algorithm estimates that relative error in the sum of squares is at most TOL.
        INFO = 2 Algorithm estimates that relatives error between X and the solution is at most TOL.
        INFO = 3 Conditions for INFO = 1 nad INFO = 2 both hold.
        INFO = 4 FVEC is orthogonal to the Jacobian to machine precision.
        INFO = 5 Number of calls to FUNC has reached or exceeded 200*(N+1).
        INFO = 6 TOL is too small, no further reduction in sum of squares is possible.
        INFO = 7 TOL is too small, no further reduction in approximate solution X is possible.
    */
    template<typename F>
    int levmarq(F&& func, std::vector<double>& x, std::vector<double>& f, double tol)
    {
        return levmarq(func, x.size(), x.data(), f.size(), f.data(), tol);
    }

}