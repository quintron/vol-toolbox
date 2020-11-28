#include <gtest/gtest.h>

#include <voltlbx/black_scholes.h>
#include <voltlbx/utils.h>
#include <cmath>

using namespace voltlbx;

TEST(BlackScholes, CheckIntrinsic)
{
    const double fwd = 3500.0;
    for(auto k: {100.0, 1000.0, 4000.0, 5000.0 })
    {        
        double zero_vol_call = bs_option_price(fwd, k, 0.0, 0.5,  1.0);
        double call_intrinsic = std::max(0.0, fwd - k); 
        ASSERT_DOUBLE_EQ(zero_vol_call, call_intrinsic);

        double zero_vol_put = bs_option_price(fwd, k, 0.0, 0.5,  -1.0);
        double put_intrinsic = std::max(0.0, k - fwd); 
        ASSERT_DOUBLE_EQ(zero_vol_put, put_intrinsic);
    }
}



TEST(BlackScholes, CheckImpliedVol)
{
    const double eps = 16.0 * std::numeric_limits<double>::epsilon();
    const double fwd = 3500.0;
    for(auto t: {1.0 / 365.0, 1.0 / 12.0, 1.0, 2.0, 5.0}) 
    {
        for(auto vol: {0.01, 0.15, 0.30, 0.50, 1.0})
        {
            const double std_dev = vol * std::sqrt(t);
            const auto nb_devs = linspace(-20.0, 20.0, 100);
            for(auto z: nb_devs)
            {
                const double strike = fwd * std::exp(z * std_dev);
                const double opt_type = z < 0.0 ? -1.0 : 1.0;
                const double price = bs_option_price(fwd, strike, vol, t, opt_type);
                const double imp_vol = bs_implied_volatility(fwd, strike, price, t, opt_type);                
                ASSERT_NEAR(imp_vol, vol, vol * eps);
            }
        }
    }
}