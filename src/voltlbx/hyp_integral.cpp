#include "hyp_integral.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <limits>
#include "utils.h"
#include "roots.h"

namespace voltlbx
{
    
    std::function<double(double)> hyperbolic_integral(double a, double b, double c)
    {
        if (c <= 0.0)
            throw std::out_of_range("c should be strictly positive");

        const double eps = std::numeric_limits<double>::epsilon();

        if (std::abs(a) <= std::abs(b) * eps) // a is null
        {
            if (std::abs(b) <= std::abs(c) * eps) // c is null
            {
                return [=](double x)-> double {return c; };
            }
            else
            {
                return [=](double x) -> double
                {
                    return 2.0 * (std::sqrt(b * x + c)- std::sqrt(c)) / b;
                };
            }
        }

        const double C = (4.0 * a * c - b * b) / (4.0 * a);
        const double sqrt_a = std::sqrt(std::abs(a));

        if (std::abs(C) == 0.0)
        {
            const double beta = 0.5 * b / a;
            return [=](double x)
            {
                return (std::log(x + beta) - std::log(beta)) / sqrt_a;
            };
        }

        const double alpha = std::sqrt(std::abs(a) / std::abs(C));
        const double beta = alpha * b * 0.5 / a;

        if (a < 0.0)
        {
            return [=](double x)
            {
                const double z = alpha * x + beta;
                return (std::asin(z) - std::asin(beta)) / sqrt_a;
            };
        }

        if (C > 0.0)
        {
            return [=](double x)
            {
                const double z = alpha * x + beta;
                return (std::asinh(z) - std::asinh(beta)) / sqrt_a;
            };
        }

        assert(C < 0.0);
        return [=](double x)
        {
            const double z = alpha * x + beta;
            if (z >= 1.0)
            {
                return (std::acosh(z) - std::acosh(beta)) / sqrt_a;
            }
            else if (z <= 1.0)
            {
                return -(std::acosh(-z) - std::acosh(-beta)) / sqrt_a;
            }
            else
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        };
    }


    std::pair<double, double> hyperbolic_integral_domain(double a, double b, double c)
    {
        if (c <= 0.0)
            throw std::out_of_range("c should be strictly positive");

        const double eps = std::numeric_limits<double>::epsilon();

        if (std::abs(a) <= std::abs(b) * eps) // a is null
        {
            if (std::abs(b) <= std::abs(c) * eps) // c is null
            {
                return { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
            }
            else
            {
                if (b > 0.0)
                    return { -c / b,  std::numeric_limits<double>::infinity() };
                else
                    return { -std::numeric_limits<double>::infinity(),  -c / b };
            }
        }

        const double C = (4.0 * a * c - b * b) / (4.0 * a);
       
        if (std::abs(C) == 0.0)
        {
            const double beta = 0.5 * b / a;
            return { -beta, std::numeric_limits<double>::infinity() };
        }

        const double alpha = std::sqrt(std::abs(a) / std::abs(C));
        const double beta = alpha * b * 0.5 / a;

        if (a < 0.0)
        {
            return { (-1.0 - beta) / alpha, (1.0 - beta) / alpha };
        }

        if (C > 0.0)
        {
            return { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        }

        assert(C < 0.0);          
        if (beta > 0.0)
        {
            return { (1.0 - beta) / alpha , std::numeric_limits<double>::infinity() };
        }

        assert(beta < 0.0);
        return { -std::numeric_limits<double>::infinity(), (-1.0 - beta) / alpha };
    }



    HyperbolicIntegral::HyperbolicIntegral(double a, double b, double c)
        :a(a),
        b(b),
        c(c),
        eval(hyperbolic_integral(a,b,c)),
        domain_(hyperbolic_integral_domain(a,b,c))
    {
        assert(c > 0.0);
    }
    

    double min_a(double b, double c, double x)
    {
        if (std::abs(x) == 0.0)
            return -std::numeric_limits<double>::infinity();

        return  -(c + x * b) / (x * x);
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