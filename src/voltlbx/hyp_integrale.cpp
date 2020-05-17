#include "hyp_integrale.h"

#include <cassert>
#include <stdexcept>
#include <limits>
#include "utils.h"
#include "roots.h"

namespace voltlbx
{       
    
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
        const double z0 = b / sqrt_disc;

        if (disc_sign < 0.0)
        {
            return (std::asinh(z) - std::asinh(z0)) / std::sqrt(a);
        }
        else if (z >= 1.0)
        {
            return (std::acosh(z) - std::acosh(z0)) / std::sqrt(a);
        }
        else if (z <= -1.0)
        {
            return (-std::acosh(-z) + std::acosh(-z0)) / std::sqrt(a);
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
                         (-1.0 - z0) * sqrt_disc / (2.0 * a) };
            }

            throw std::logic_error("should never get there !");
        }

        throw std::runtime_error("not yet implemented");
    }


    double min_a(double b, double c, double x)
    {
        const double min_val = c * 252.0 * std::numeric_limits<double>::epsilon();
        
        if (std::abs(x) == 0.0)
            return min_val;

        return std::max(min_val, -(c + x * b) / (x * x));
    }


    double HyperbolicIntegral::max_value(double b, double c, double x)
    {            
        return HyperbolicIntegral(min_a(b, c, x), b, c)(x);
    }


    HyperbolicIntegral HyperbolicIntegral::solve_a_from_value(double b, double c, double x, double value)
    {
        if (c <= 0.0)
            throw std::out_of_range("c should be strictly positive");           
        
        if(x * value < 0.0)
            throw std::out_of_range("value and x should have same sign");
        
        if ( std::abs(max_value(b, c, x)) < std::abs(value))
            throw std::runtime_error("value is not attainable");

        const double min_search_a = min_a(b, c, x);
        assert(min_search_a > 0.0);
        
        const auto err = [&](double u)
        {
            double a = min_search_a + std::exp(u);
            return HyperbolicIntegral(a, b, c)(x) - value;
        };
        const auto[u1, u2] = bracket_root(err, -10.0, 0.0);
        const double u = brenth(err, u1, u2);
        const double a = min_search_a + std::exp(u);
        return HyperbolicIntegral(a, b, c);
    }

}