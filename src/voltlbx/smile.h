#pragma once
#include "pimpl.h"
#include "jet.h"
#include <cmath>
#include <vector>

namespace voltlbx
{
    class Smile 
    {
    public:

        Smile(double time_to_maturity, double atf_vol)
            :time_to_maturity(time_to_maturity), 
            atf_vol(atf_vol),             
            atf_dev(atf_vol * std::sqrt(this->time_to_maturity))
        {            
        }

        const double time_to_maturity;
        const double atf_vol;                

        virtual ~Smile(){}
        
        virtual double backbone(double z) const = 0;

        virtual Jet backbone_jet(double z) const = 0;

        virtual double vol(double x) const
        {
            return atf_vol * backbone(x / atf_dev);
        }

        /// <summary>
        /// Ratio of smile density and black-scholes density. 
        /// </summary>
        double density_ratio(double x) const;


        /// <summary>
        /// BBF (Beresticky-Busca-Florent) pseudo local vol : smile volatility is harmonic average of pseudo lv        /// 
        /// </summary>
        double pseudo_local_vol(double x) const;

    private:
        const double atf_dev;
    };

    
    class CubicSplineSmile :  public Smile, Pimpl<CubicSplineSmile>
    {
    public :
        CubicSplineSmile(double time_to_maturity, 
                         double atf_vol,
                         std::vector<double> zs,
                         std::vector<double> vol_ratios);

        double backbone(double z) const override;

        Jet backbone_jet(double z) const override;
    };

}