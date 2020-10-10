#include "smile.h"

#include "cubic_spline.h"
#include "utils.h"
#include <algorithm>
#include <cassert>

namespace voltlbx
{
    double density_ratio_from_backbone(double z, double atf_dev, double b, double db_dz, double d2b_d2z)
    {
        return (1.0 - z * db_dz / b + 0.5 * atf_dev * b * db_dz)
                * (1.0 - z * db_dz / b - 0.5 * atf_dev * b * db_dz)
               + b * d2b_d2z;
    }


    double Smile::density_ratio(double x) const
    {        
        const double z = x / atf_dev;
        auto j = backbone_jet(z);
        return density_ratio_from_backbone(z, atf_dev, j.y, j.dy_dx, j.d2y_d2x);  
    }


    double pseudo_lv(double z, double atf_vol, double b, double db_dz)
    {
        return atf_vol * b / (1.0 - z * db_dz / b);
    }


    double Smile::pseudo_local_vol(double x) const
    {
        const double z = x / atf_dev;
        auto j = backbone_jet(z);
        return pseudo_lv(z, atf_vol, j.y, j.dy_dx);
    }


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
        const Jet left_sqr_jet;
        
        const double right_z;
        const Jet right_sqr_jet;
    };

    CubicSplineSmile::CubicSplineSmile(double time_to_maturity,
                                       double atf_vol, 
                                       std::vector<double> zs, 
                                       std::vector<double> vol_ratios)
        : Smile(time_to_maturity, atf_vol), 
        Pimpl<CubicSplineSmile>(zs, vol_ratios)
    {
    }


    Jet CubicSplineSmile::backbone_jet(double z) const
    {
        const auto sqr_ratio = [&]()
        {
            if (z < impl->left_z)
            {
                const auto [y, yp, _] = impl->left_sqr_jet;
                return Jet{ y + yp * (z - impl->left_z), yp, 0.0 };
            }

            if (z > impl->right_z)
            {
                const auto [y, yp, _] = impl->right_sqr_jet;
                return Jet{ y + yp * (z - impl->right_z), yp, 0.0 };
            }

            return impl->sqr_ratio.eval_jet(z);
        }();

        return pow(sqr_ratio, 0.5);
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
                const auto[y, yp, _] = impl->right_sqr_jet;
                return y + yp * (z - impl->right_z);
            }

            return impl->sqr_ratio(z);
        }();

        return std::sqrt(std::max(0.0, sqr_ratio));
    }
}