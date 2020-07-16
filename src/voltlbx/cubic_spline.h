#pragma once

#include <vector>
#include <optional>

namespace voltlbx
{
    /// Cubic spline interpolation
    class CubicSpline
    {
    public:

        /// <summary>
        /// Build cubic spline. Algorithm from numerical recipes.
        /// 
        /// Left boundary condition  : 
        ///     spline derivative match left_derivative when it has value
        ///     otherwise "natural" spline condition (zero second order derivative)        
        /// 
        /// Right boundary condition  : 
        ///     spline derivative match right_derivative when it has value
        ///     otherwise "natural" spline condition (zero second order derivative)   
        /// </summary>
        /// <param name="xs"> abscissae nodes (assumed to be sorted) </param>
        /// <param name="ys"> values nodes </param>
        /// <param name="left_derivative">Left boundary condition, default is "natural" : zero second derivative </param>
        /// <param name="right_derivative">Right boundary condition, default is "natural" : zero second derivative </param>
        CubicSpline(std::vector<double> xs,
                    std::vector<double> ys,
                    std::optional<double> left_derivative,
                    std::optional<double> right_derivative);

        /// <summary>
        /// Build natural cubic spline. 
        /// </summary>
        /// <param name="xs"> abscissae nodes (assumed to be sorted) </param>
        /// <param name="ys"> values nodes </param>
        CubicSpline(std::vector<double> xs,
                    std::vector<double> ys)
            :CubicSpline(std::move(xs), std::move(ys), {}, {})
        {
        }

        double operator()(double) const;

        double eval(double x, double* yp, double* ypp) const;
        std::tuple<double, double, double> eval_jet(double x) const
        {
            double yp, ypp;
            double y = eval(x, &yp, &ypp);
            return { y, yp, ypp };
        }

    private:
        const std::vector<double> xs_;
        const std::vector<double> ys_;
        std::vector<double> ypps;
    };


}