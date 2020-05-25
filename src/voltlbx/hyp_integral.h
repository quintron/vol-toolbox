#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include <map>

namespace voltlbx
{

    /* Primitive of 1 / sqrt(a * x^2 + b * x + c)
       Branch that contains zero and is vanishing at 0.
       Contraint is that c > 0 and a > 0.
    */
    class HyperbolicIntegral
    {
    public:
        HyperbolicIntegral(double a, double b, double c);

        static HyperbolicIntegral solve_a_from_value(double b, double c, double x, double value);
        static double max_value(double b, double c, double x);

        double operator()(double x) const;

        std::pair<double, double> domain() const;
        
        const double a;
        const double b;
        const double c;

    private:
        const double sqrt_disc;
        const int disc_sign;
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