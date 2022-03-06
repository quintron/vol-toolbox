#pragma once
#include <vector>

#include "pimpl.h"

namespace voltlbx
{
    class SmileVariationFilter : public Pimpl<SmileVariationFilter>
    {
    public:
        struct Config
        {
            double atm_dev;
            double atm_skew_ratio;
            double wing_skew_ratio;
            double z_ref;
        };

        SmileVariationFilter(const std::vector<double>& zs,
                             const std::vector<double>& dvols,
                             const std::vector<double>& error_devs, 
                             Config config);

        static SmileVariationFilter create(const std::vector<double>& zs,
                                           const std::vector<double>& dvols,
                                           const std::vector<double>& error_devs, 
                                           double atm_dev,
                                           double atm_skew_ratio,
                                           double wing_skew_ratio,
                                           double z_ref)
        {
            return SmileVariationFilter(zs, dvols, error_devs, 
                                        { atm_dev, atm_skew_ratio, wing_skew_ratio, z_ref });
        }

        std::tuple<double, double> dvol_with_error(double z) const;

        double dvol(double z) const;

    };
}
