#include "hyp_integral.h"

#include <cassert>
#include <cmath>
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
        const auto[u1, u2] = bracket_root(err, -2.0, 0.0);
        const double u = brenth(err, u1, u2);
        const double a = min_search_a + std::exp(u);
        return HyperbolicIntegral(a, b, c);
    }


    double MonotoneHyperbolicSpline::operator()(double x) const
    {        
        const int i = std::max(0, std::min(static_cast<int>(xs.size()) - 2, locate_left_index(xs, x)));        
        return step_funcs[i](x - xs[i]) + ys[i];
    }

    MonotoneHyperbolicSpline MonotoneHyperbolicSpline::build
        (std::map<double, double> nodes,
         double left_slope,
         double left_convexity)
    {
        std::vector<double> xs;
        std::vector<double> ys;
        for (auto[x, y] : nodes)
        {
            xs.push_back(x);
            ys.push_back(y);
        }

        // 1 / sqrt(ax^2 + bx + c )        
        // -(ax + b / 2) / (ax^2 + bx + c)^(3/2)

        double current_y = ys[0];
        double current_c = 1.0 / (left_slope * left_slope);
        double current_b = -2.0 * left_convexity * std::pow(current_c, 1.5);        
        std::vector<HyperbolicIntegral> step_funcs;
        for (std::size_t i = 0; i + 1 < xs.size(); ++i)
        {
            const double d_y = ys[i + 1] - current_y;
            if (d_y <= 0.0)
            {
                throw std::logic_error("values should be monotone");
            }

            const double d_x = xs[i + 1] - xs[i];
            const auto step_f = HyperbolicIntegral::solve_a_from_value(current_b, current_c, d_x, d_y);
            step_funcs.push_back(step_f);

            current_y += step_f(d_x);
            current_c = step_f.a * d_x * d_x + step_f.b * d_x + step_f.c;
            current_b = 2.0 * step_f.a * d_x + step_f.b;
        }

        return MonotoneHyperbolicSpline(xs, ys, step_funcs);
    }

}