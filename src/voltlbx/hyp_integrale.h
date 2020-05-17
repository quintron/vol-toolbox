#pragma once

#include <cmath>
#include <utility>

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
    
}