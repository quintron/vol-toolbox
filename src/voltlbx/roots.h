#pragma once

#include <functional>

namespace voltlbx
{           
    std::pair<double, double> bracket_root(const std::function<double(double)>& f, double x1, double x2);

   
    /// <summary>
    /// Brenth algorithm from scipy
    /// </summary>
    double brenth(const std::function<double(double)>& f, 
                  double xa, double xb,
                  double xtol = 2e-12, double rtol = 8.881784197001252e-16, std::size_t max_iter = 100);
}