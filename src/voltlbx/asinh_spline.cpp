#include "asinh_spline.h"

#include <cassert>
#include <stdexcept>
#include "utils.h"
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


    HyperbolicIntegral::HyperbolicIntegral(double a, double b, double c)
        :a(a),
        b(b),
        c(c),
        sqrt_disc(std::sqrt(std::abs(b * b - 4.0 * a * c))),
        disc_sign(sgn(b*b - 4.0 * a * c))
    {
        assert(a > 0.0);
        assert(c > 0.0);
    }

    double HyperbolicIntegral::operator()(double x) const
    {
        if (disc_sign == 0)
        {
            const double u = x + b / (2.0 * a);
            return std::log(u) / std::sqrt(a);
        }

        const double z = (2.0 * a * x + b) / sqrt_disc;

        if (disc_sign < 0.0)
        {
            return std::asinh(z) / std::sqrt(a);
        }
        else if (z >= 1.0)
        {
            return std::acosh(z) / std::sqrt(a);
        }
        else if (z <= -1.0)
        {
            return -std::acosh(-z) / std::sqrt(a);
        }
        else
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::pair<double, double> HyperbolicIntegral::domain() const
    {
        if (disc_sign < 0.0)
        {
            return { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        }
        else if (disc_sign > 0.0)
        {
            double z0 = b / sqrt_disc;
            if (z0 >= 1.0)
            {
                return { (1.0 - z0) * sqrt_disc / (2.0 * a),
                         std::numeric_limits<double>::infinity() };
            }
            else if (z0 <= -1.0)
            {
                return { -std::numeric_limits<double>::infinity(),
                         (-1.0 - z0) *sqrt_disc / (2.0 * a) };
            }

            throw std::logic_error("should never get there !");
        }

        throw std::runtime_error("not yet implemented");
    }

}