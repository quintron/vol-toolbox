#pragma once

#include <functional>

namespace voltlbx
{
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

    template<typename F>
    int levmarq(F&& func, std::vector<double>& x, std::vector<double>& f, double tol)
    {
        return levmarq(func, x.size(), x.data(), f.size(), f.data(), tol);
    }

}