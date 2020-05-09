#pragma once

#include <cmath>

namespace voltlbx
{
     
    class AsinhSplineElem
    {
    public:
        AsinhSplineElem(double alpha, double beta)
            :alpha(alpha),
            beta(beta),
            c1(std::sqrt(4.0 * alpha - beta * beta)),
            c2(std::sqrt(alpha))
        {
        }

        double operator()(double u) const
        {            
            return (std::asinh((2.0 * alpha * u + beta) / c1) - std::asinh(beta / c1)) / c2;            
        }

    private:
        const double alpha;
        const double beta;
        const double c1;
        const double c2;
    };

    double max_value1_from_left_second_deriv(double second_deriv0);

    AsinhSplineElem solve_from_left_second_deriv_and_right_value(double second_deriv0, double  value1);
    
}