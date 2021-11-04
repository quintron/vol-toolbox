#include "smile.h"

#include "cubic_spline.h"
#include "utils.h"
#include <algorithm>
#include <cassert>

namespace voltlbx
{

    template<>
    struct Pimpl<ExtrapolatedSplineCurve>::Implementation
    {
        Implementation(std::vector<double> xs,
                       std::vector<double> values)
            : interpol(xs, values),
            left_x(xs.front()),
            left_jet(interpol.eval_jet(left_x)),
            right_x(xs.back()),
            right_jet(interpol.eval_jet(right_x))
        {
        }

        const CubicSpline interpol;

        const double left_x;
        const Jet left_jet;

        const double right_x;
        const Jet right_jet;
    };

    ExtrapolatedSplineCurve::ExtrapolatedSplineCurve(std::vector<double> xs,
                                                     std::vector<double> values)
        : Pimpl<ExtrapolatedSplineCurve>(xs, values)
    {
    }

    Jet ExtrapolatedSplineCurve::vol_jet(double x) const
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


    template<>
    struct Pimpl<SplineSmileCurve>::Implementation
    {
        Implementation(std::vector<double> xs,
                       std::vector<double> vols)
            : sqr_vol(xs, map(vols, [](double v) {return v * v; }))
        {
        }

        const ExtrapolatedSplineCurve sqr_vol;
    };


    SplineSmileCurve::SplineSmileCurve(std::vector<double> xs,
                                       std::vector<double> vols)
        : Pimpl<SplineSmileCurve>(xs, vols)
    {
    }


    Jet SplineSmileCurve::vol_jet(double x) const
    {
        const Jet sqr_v = impl->sqr_vol.vol_jet(x);
        return pow(sqr_v, 0.5);
    }

}