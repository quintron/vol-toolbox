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


    class SplineSmileCurve : public SmileCurve, Pimpl<SplineSmileCurve>
    {
    public:
        SplineSmileCurve(std::vector<double> xs,
                         std::vector<double> vols);

        Jet vol_jet(double x) const override;
    };

}