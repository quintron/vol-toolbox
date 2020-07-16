#include "smile.h"

#include "cubic_spline.h"
#include "utils.h"
#include <algorithm>

namespace voltlbx
{
    template<>
    struct Pimpl<CubicSplineSmile>::Implementation
    {    
        Implementation(std::vector<double> zs,
                       std::vector<double> vol_ratios)
            : sqr_ratio(zs, map(vol_ratios, [](double r) {return r * r; })),
            left_z(zs.front()),            
            left_sqr_jet(sqr_ratio.eval_jet(left_z)),
            right_z(zs.back()),
            right_sqr_jet(sqr_ratio.eval_jet(right_z))
        {
        }

        const CubicSpline sqr_ratio;

        const double left_z;
        const std::tuple<double, double, double> left_sqr_jet;
        
        const double right_z;
        const std::tuple<double, double, double> right_sqr_jet;
    };

    CubicSplineSmile::CubicSplineSmile(double time_to_maturity,
                                       double atf_vol, 
                                       std::vector<double> zs, 
                                       std::vector<double> vol_ratios)
        : Smile(time_to_maturity, atf_vol), 
        Pimpl<CubicSplineSmile>(zs, vol_ratios)
    {
    }

    double CubicSplineSmile::backbone(double z) const
    {        
        const double sqr_ratio = [&]()
        {
            
            if (z < impl->left_z)
            {
                const auto[y, yp, _] = impl->left_sqr_jet;
                return y + yp * (z - impl->left_z);
            }

            if (z > impl->right_z)
            {
                const auto [y, yp, _] = impl->right_sqr_jet;
                return y + yp * (z - impl->right_z);
            }

            return impl->sqr_ratio(z);
        }();

        return std::sqrt(std::max(0.0, sqr_ratio));
    }
}