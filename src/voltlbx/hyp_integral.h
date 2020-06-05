#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include <map>
#include <functional>

namespace voltlbx
{

    /* Primitive of 1 / sqrt(a * x^2 + b * x + c)
       Branch defined and vanishing at 0 (it is then required that c > 0).
    */
    class HyperbolicIntegral
    {
    public:
        HyperbolicIntegral(double a, double b, double c);

        static HyperbolicIntegral solve_a_from_value(double b, double c, double x, double value);
        static double max_value(double b, double c, double x);
                
        double operator()(double x) const
        {
            return eval(x);
        }

        std::pair<double, double> domain() const
        {
            return domain_;
        }
        
        const double a;
        const double b;
        const double c;

    private:
        const std::function<double(double)> eval;
        const std::pair<double, double> domain_;
    };
    


    class MonotoneHyperbolicSpline
    {
    public:

        static MonotoneHyperbolicSpline build
            (std::map<double, double> nodes,
             double left_slope,
             double left_convexity);

        double operator()(double z) const;

    private:
        MonotoneHyperbolicSpline(std::vector<double> xs,
                                 std::vector<double> ys,
                                 std::vector<HyperbolicIntegral> step_funcs)
            :xs(std::move(xs)),
            ys(std::move(ys)),
            step_funcs(std::move(step_funcs))
        {
        }

        const std::vector<double> xs;
        const std::vector<double> ys;
        const std::vector<HyperbolicIntegral> step_funcs;
    };
          
}