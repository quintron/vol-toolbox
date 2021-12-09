#pragma once
#include "pimpl.h"
#include "jet.h"
#include <cmath>
#include <vector>

namespace voltlbx
{
    
    class SmileCurve
    {
    public:

        virtual Jet vol_jet(double x) const = 0;
        
        double vol(double x) const
        {
            return vol_jet(x).y;
        }
    };


    class SplineCurve : public SmileCurve, Pimpl<SplineCurve>
    {
    public:
        SplineCurve(const std::vector<double>& xs,
                    const std::vector<double>& values);

        Jet vol_jet(double x) const override;

        const std::vector<double>& nodes() const;
    };


    class SmileSplineCurve : public SmileCurve, Pimpl<SmileSplineCurve>
    {
    public:
        SmileSplineCurve(std::vector<double> xs,
                         std::vector<double> vols);

        Jet vol_jet(double x) const override;
    };


    class AverageSplineCurve : public SmileCurve, Pimpl<AverageSplineCurve>
    {
    public:
        AverageSplineCurve(std::vector<double> xs,
                           std::vector<double> values);

        Jet vol_jet(double x) const override;
    };

}
