#include "smile.h"

#include "cubic_spline.h"
#include "utils.h"
#include "quad.h"
#include <algorithm>
#include <cassert>

namespace voltlbx
{

    template<>
    struct Pimpl<SplineCurve>::Implementation
    {
        Implementation(const std::vector<double>& xs,
                       const std::vector<double>& values)
            :
            nodes(xs),
            interpol(xs, values),
            left_x(xs.front()),
            left_jet(interpol.eval_jet(left_x)),
            right_x(xs.back()),
            right_jet(interpol.eval_jet(right_x))
        {
        }

        const std::vector<double> nodes;
        const CubicSpline interpol;

        const double left_x;
        const Jet left_jet;

        const double right_x;
        const Jet right_jet;
    };

    SplineCurve::SplineCurve(const std::vector<double>& xs,
                                         const std::vector<double>& values)
        : Pimpl<SplineCurve>(xs, values)
    {
    }


    Jet SplineCurve::vol_jet(double x) const
    {
        if (x < impl->left_x)
        {
                const auto [y, yp, _] = impl->left_jet;
                return Jet{ y + yp * (x - impl->left_x), yp, 0.0 };
        }

        if (x > impl->right_x)
        {
                const auto [y, yp, _] = impl->right_jet;
                return Jet{ y + yp * (x - impl->right_x), yp, 0.0 };
        }

        return impl->interpol.eval_jet(x);        
    }


    const std::vector<double>& SplineCurve::nodes() const
    {
        return impl->nodes;
    }


    template<>
    struct Pimpl<SmileSplineCurve>::Implementation
    {
        Implementation(std::vector<double> xs,
                       std::vector<double> vols)
            : sqr_vol(xs, map(vols, [](double v) {return v * v; }))
        {
        }

        const SplineCurve sqr_vol;
    };


    SmileSplineCurve::SmileSplineCurve(std::vector<double> xs,
                                       std::vector<double> vols)
        : Pimpl<SmileSplineCurve>(xs, vols)
    {
    }


    Jet SmileSplineCurve::vol_jet(double x) const
    {
        const Jet sqr_v = impl->sqr_vol.vol_jet(x);
        return pow(sqr_v, 0.5);
    }


    double quad_integral(double a, double b, const SplineCurve& crv)
    {
        auto quad = GaussLegendreInterval<4>(a, b);
        return quad([&](double x) { return crv.vol(x); });
    }

    template<>
    struct Pimpl<AverageSplineCurve>::Implementation
    {
        Implementation(std::vector<double> xs,
                       std::vector<double> values)
            : inner_crv(xs, values)            
        {
            const auto& nodes = inner_crv.nodes();
            i0 = static_cast<size_t>(
                 std::min(static_cast<int>(nodes.size()) - 1,
                          std::max(0, locate_left_index(nodes, 0.0))));
            int0 = quad_integral(0.0, nodes[i0], inner_crv);
            for (std::size_t i = 0; i < nodes.size(); i++)
            {
                ints.push_back(quad_integral(nodes[i], nodes[i+1], inner_crv));
            }
        }
         
        double integral(double z) const
        {
            const auto& nodes = inner_crv.nodes();
            auto iz = static_cast<size_t>(
                        std::min(static_cast<int>(nodes.size()) - 1,
                            std::max(0, locate_left_index(nodes, z))));

            double res = int0;
            res += quad_integral(nodes[iz], z, inner_crv);

            if (i0 < iz)
            {   
                for (std::size_t i = i0; i < iz; ++i)
                {
                    res += ints[i];
                }
            }
            else if(iz < i0)
            {
                for (std::size_t i = iz; i < i0; ++i)
                {
                    res -= ints[i];
                }
            }

            return res;
        }

        double average(double z) const
        {
            constexpr double eps = std::numeric_limits<double>::epsilon();
            if (std::abs(z) <= eps)
            {
                return inner_crv.vol(0.0);
            }

            return integral(z) / z;
        }

        const SplineCurve inner_crv;

        size_t i0;
        double int0;
        std::vector<double> ints;

    };


    AverageSplineCurve::AverageSplineCurve(std::vector<double> xs,
        std::vector<double> values)
        : Pimpl<AverageSplineCurve>(xs, values)
    {
    }


    Jet AverageSplineCurve::vol_jet(double x) const
    {
        return Jet 
        { impl->average(x), 
            std::numeric_limits<double>::quiet_NaN(), // TODO
          std::numeric_limits<double>::quiet_NaN() }; // TODO
    }

}
