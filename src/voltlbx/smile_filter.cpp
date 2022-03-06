#include "smile_filter.h"

#include <cmath>
#include "linear_filter.h"

namespace voltlbx
{
    class SmileCov : public Covariance<double>
    {
    public:

        SmileCov(double atm_dev,
                 double atm_skew_ratio,
                 double wing_skew_ratio,
                 double z_ref)
            :atm_dev(atm_dev),
            atm_skew_ratio(atm_skew_ratio),
            wing_skew_ratio(wing_skew_ratio),
            z_ref(z_ref)
        {
        }

        double dev(double z) const override
        {
            return atm_dev * std::sqrt(1.0 + std::pow(wing_skew_ratio * z, 2));
        }

        double correl(double z1, double z2) const override
        {
            return std::exp(-0.5 * std::pow(atm_skew_ratio * z_ref * (std::atan(z1 / z_ref) - std::atan(z2 / z_ref)), 2.0));
        }

    private:
        const double atm_dev;
        const double atm_skew_ratio;
        const double wing_skew_ratio;
        const double z_ref;
    };


    template<>
    struct Pimpl<SmileVariationFilter>::Implementation
    {
        explicit Implementation(const std::vector<double>& zs,
                                const std::vector<double>& dvols,
                                const std::vector<double>& error_devs,
                                double atm_dev,
                                double atm_skew_ratio,
                                double wing_skew_ratio,
                                double z_ref)
            : filter(std::make_shared<SmileCov>(atm_dev, atm_skew_ratio, wing_skew_ratio, z_ref))
        {
            filter.fit(zs, dvols, error_devs);
        }

        std::tuple<double, double> dvol_with_error(double z) const
        {
            return filter.predict_with_error(z);
        }

        double dvol(double z) const
        {
            return filter.predict(z);
        }

        LinearFilter<double> filter;
    };


    SmileVariationFilter::SmileVariationFilter(const std::vector<double>& zs,
                                               const std::vector<double>& dvols,
                                               const std::vector<double>& error_devs, 
                                               Config config)
        : Pimpl<SmileVariationFilter>(zs, dvols, error_devs,
                                      config.atm_dev,
                                      config.atm_skew_ratio,
                                      config.wing_skew_ratio,
                                      config.z_ref)
    {
    }


    std::tuple<double, double> SmileVariationFilter::dvol_with_error(double z) const
    {
        return impl->dvol_with_error(z);
    }


    double SmileVariationFilter::dvol(double z) const
    {
        return impl->dvol(z);
    }

}
