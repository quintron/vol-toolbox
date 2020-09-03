#pragma once
#include <cmath>
#include <cassert>

namespace voltlbx
{
    struct Jet
    {
        double y;
        double dy_dx;
        double d2y_d2x;

        double operator()(double dy) const
        {
            return y + (dy_dx + 0.5 * d2y_d2x * dy) * dy;
        }
    };


    inline Jet operator+(const Jet& f, const Jet& g)
    {
        return Jet{ f.y + g.y, f.dy_dx + g.dy_dx, f.d2y_d2x + g.d2y_d2x };
    }

    inline Jet operator+(const Jet& j, const double& a)
    {
        return Jet{ a + j.y, j.dy_dx, j.d2y_d2x };
    }

    inline Jet operator+(const double& a, const Jet& j)
    {
        return j + a;
    }

    inline Jet operator*(const Jet& f, const Jet& g)
    {
        return Jet{ f.y * g.y, f.dy_dx * g.y + f.y * g.dy_dx, f.d2y_d2x * g.y + g.d2y_d2x * f.y + 2.0 * f.dy_dx * g.dy_dx };
    }

    inline Jet operator*(const Jet& j, const double& a)
    {
        return Jet{ a * j.y, a * j.dy_dx, a * j.d2y_d2x };
    }

    inline Jet operator*(const double& a, const Jet& j)
    {
        return j * a;
    }


    inline Jet pow(const Jet& j, double p)
    {
        assert(j.y > 0.0);
        //(1 + u)^p = 1 + p *u + 0.5 * p * (p - 1 ) * u^2
        auto norm_j = Jet{ 0.0, j.dy_dx / j.y,  j.d2y_d2x / j.y };
        auto norm_j_p = 1.0 + p * norm_j + (0.5 * p * (p - 1.0)) * (norm_j * norm_j);
        return std::pow(j.y, p) * norm_j_p;
    }

}