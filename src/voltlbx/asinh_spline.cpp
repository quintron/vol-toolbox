#include "asinh_spline.h"

#include <cassert>
#include <stdexcept>

#include "roots.h"

namespace voltlbx
{       
    double max_value1_from_left_second_deriv(double second_deriv0)
    {
        assert(second_deriv0 < 4.0);
        const double beta = -0.5 * second_deriv0;
        return 2.0 * std::log(beta / 2.0 + 1.0) / beta;
    }

    AsinhSplineElem solve_from_left_second_deriv_and_right_value(double second_deriv0, double value1)
    {
        assert(second_deriv0 < 4.0);
        const double beta = -0.5 * second_deriv0;        
        const double min_alpha = 0.25 * beta * beta;

        const auto err = [&](double x)
        {
               double alpha = min_alpha + std::exp(x);
               return AsinhSplineElem(alpha, beta)(1.0) - value1;
        };
        const auto[x1, x2] = bracket_root(err, -10.0, std::log(min_alpha));
        const double x = brenth(err, x1, x2);
        const double alpha = min_alpha + std::exp(x);
        return AsinhSplineElem(alpha, beta);
    }    

}