#include "smile.h"

#include "cubic_spline.h"
#include "utils.h"
#include <algorithm>
#include <cassert>

namespace voltlbx
{

    template<>
    struct Pimpl<SplineSmileCurve>::Implementation
    {
        Implementation(std::vector<double> xs,
                       std::vector<double> vols)
            : sqr_ratio(xs, map(vols, [](double v) {return v * v; })),
            left_x(xs.front()),
            left_sqr_jet(sqr_ratio.eval_jet(left_x)),
            right_x(xs.back()),
            right_sqr_jet(sqr_ratio.eval_jet(right_x))
        {
        }

        const CubicSpline sqr_ratio;

        const double left_x;
        const Jet left_sqr_jet;

        const double right_x;
        const Jet right_sqr_jet;
    };


    SplineSmileCurve::SplineSmileCurve(std::vector<double> xs,
                                       std::vector<double> vols)
        : Pimpl<SplineSmileCurve>(xs, vols)
    {
    }


    Jet SplineSmileCurve::vol_jet(double x) const
    {
        const auto sqr_ratio = [&]()
        {
            if (x < impl->left_x)
            {
                const auto [y, yp, _] = impl->left_sqr_jet;
                return Jet{ y + yp * (x - impl->left_x), yp, 0.0 };
            }

            if (x > impl->right_x)
            {
                const auto [y, yp, _] = impl->right_sqr_jet;
                return Jet{ y + yp * (x - impl->right_x), yp, 0.0 };
            }

            return impl->sqr_ratio.eval_jet(x);
        }();

        return pow(sqr_ratio, 0.5);
    }

}