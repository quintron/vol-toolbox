#include "cubic_spline.h"

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include "utils.h"

namespace voltlbx
{
    
    CubicSpline::CubicSpline(std::vector<double> xs, 
                             std::vector<double> ys, 
                             std::optional<double> left_derivative, 
                             std::optional<double> right_derivative)
        :xs_(std::move(xs)),
        ys_(std::move(ys)),
        ypps(std::vector<double>(xs_.size()))
    {
        if (xs_.empty())
            throw std::length_error("empty input not allowed");

        if (xs_.size() != ys_.size())
            throw std::length_error("inputs should have same size");
                
        std::vector<double> u(xs_.size());

        if (!left_derivative.has_value()) //Natural Spline
        {
            ypps[0] = 0.0;
        }
        else
        {
            ypps[0] = -0.5;
            u[0] = (3.0 / (xs_[1] - xs_[0])) * ((ys_[1] - ys_[0]) / (xs_[1] - xs_[0]) - left_derivative.value());
        }

        for (std::size_t i = 1; i < xs_.size() - 1; i++)
        {
            double sig = (xs_[i] - xs_[i - 1]) / (xs_[i + 1] - xs_[i - 1]);
            double p = sig * ypps[i - 1] + 2.0;
            ypps[i] = (sig - 1.0) / p;
            u[i] = (ys_[i + 1] - ys_[i]) / (xs_[i + 1] - xs_[i]) - (ys_[i] - ys_[i - 1]) / (xs_[i] - xs_[i - 1]);
            u[i] = (6.0 * u[i] / (xs_[i + 1] - xs_[i - 1]) - sig * u[i - 1]) / p;
        }

        int n = xs_.size();
        double qn, un;
        if (!right_derivative.has_value()) //Natural Spline
        {
            qn = 0.0;
            un = 0.0;
        }
        else
        {
            qn = 0.5;
            un = (3.0 / (xs_[n - 1] - xs_[n - 1])) *
                (right_derivative.value() - (ys_[n - 1] - ys_[n - 2]) / (xs_[n - 1] - xs_[n - 2]));
        }
        ypps[n - 1] = (un - qn * u[n - 2]) / (qn * ypps[n - 2] + 1.0);
        for (int k = n - 2; k >= 0; k--)
        {
            ypps[k] = ypps[k] * ypps[k + 1] + u[k];
        }
    }
    
    double CubicSpline::operator()(double x) const
    {
        const int i = std::min(static_cast<int>(xs_.size()) - 2, 
                               std::max(0, locate_left_index(xs_, x)));

        const double h = xs_[i + 1] - xs_[i];
        assert(h > 0.0);
        const double a = (xs_[i + 1] - x) / h;
        const double b = (x - xs_[i]) / h;        
        return a * ys_[i] + b * ys_[i + 1] + (a * (a * a - 1.0) * ypps[i] + b * (b * b - 1.0) * ypps[i + 1]) * (h * h) / 6.0;
    }

}