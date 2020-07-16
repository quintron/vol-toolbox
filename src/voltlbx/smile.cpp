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
            : sqr_ratio(zs, map(vol_ratios, [](double r) {return r * r; }))
        {
        }

        const CubicSpline sqr_ratio;
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
        return std::sqrt(std::max(0.0, impl->sqr_ratio(z)));
    }
}