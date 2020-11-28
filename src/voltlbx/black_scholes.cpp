#include "black_scholes.h"

#include "LetsBeRational/lets_be_rational.h"

namespace voltlbx
{
    double bs_option_price(double forward, double strike, double vol, double time, double option_type)
    {
        return black(forward, strike, vol, time, option_type);
    }    

    double bs_implied_volatility(double forward, double strike, double price, double time, double option_type)
    {
        return implied_volatility_from_a_transformed_rational_guess(price, forward, strike, time, option_type);
    }
}