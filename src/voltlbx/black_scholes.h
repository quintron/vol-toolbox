#pragma once

#include <vector>
#include <optional>

namespace voltlbx
{

    double bs_option_price(double forward, double strike, double vol, double time, double option_type);

    double bs_implied_volatility(double forward, double strike, double price, double time, double option_type);

}