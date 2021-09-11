#pragma once
#include "pimpl.h"
#include <Eigen/Dense>
#include <stdexcept>

namespace voltlbx
{

    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;
    using SVD = Eigen::BDCSVD<Eigen::Matrix<double, -1, -1, 0>>;
    
    template<typename T>
    class Covariance
    {
    public:

        virtual double dev(T) const = 0;

        virtual double correl(T, T) const = 0;

        virtual double operator()(T a, T b) const
        {
            return correl(a, b) * dev(a) * dev(b);
        }
    };

    template<typename T>
    class LinearFilter
    {
    public:

        using Cov = std::shared_ptr<Covariance<T>>;

        LinearFilter(Cov covariance)
            :covariance(covariance)
        {}

        void fit(
            const std::vector<T>& inputs,
            const std::vector<double>& targets,
            const std::vector<double>& error_devs)
        {
            if (inputs.size() != targets.size()
                || inputs.size() != error_devs.size())
            {
                throw std::logic_error("LinearFilter inputs should have same size");
            }

            _inputs_devs = std::vector<double>(inputs.size());            
            for (std::size_t i = 0; i < inputs.size(); ++i)
            {
                _inputs_devs[i] = (*covariance).dev(inputs[i]);
            }

            Matrix cov = Matrix::Zero(inputs.size(), inputs.size());
            for (std::size_t i = 0; i < inputs.size(); ++i)
            {
                const double di = _inputs_devs[i];                    
                for (std::size_t j = 0; j < i; ++j)
                {
                    cov(i, j) = di * _inputs_devs[j] * (*covariance).correl(inputs[i], inputs[j]);
                    cov(j, i) = cov(i, j);
                }
                cov(i, i) = di * di * (*covariance).correl(inputs[i], inputs[i]);
            }
            for (std::size_t i = 0; i < inputs.size(); ++i)
            {
                const double e = error_devs[i];
                cov(i, i) += e * e;                    
            }

            cov_svd = cov.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
            const auto v_targets = Vector::Map(targets.data(), targets.size());
            weights = cov_svd.solve(v_targets);
            _inputs = inputs;
        }

        std::pair<double, double> predict_with_error(T new_input, bool with_error = true) const
        {
            if (_inputs.empty())
                throw std::logic_error("Linear filter should be fitted before prediction");

            Vector cov = cov_with_inputs(new_input);
            const double prediction = cov.dot(weights);

            double err = std::numeric_limits<double>::quiet_NaN();
            if (with_error)
            {
                err = std::sqrt(std::max(0.0, (*covariance)(new_input, new_input) - cov.dot(cov_svd.solve(cov))));
            }
            return { prediction, err };
        }

        double predict(T new_input) const
        {
            auto[p, e] = predict_with_error(new_input, false);
            return p;
        }

    private:

        Vector cov_with_inputs(T new_input) const
        {
            const double v_input = (*covariance).dev(new_input);
            Vector cov_with_inputs = Vector::Zero(_inputs.size());
            for (std::size_t i = 0; i < _inputs.size(); ++i)
            {
                cov_with_inputs[i] = v_input * _inputs_devs[i] * (*covariance).correl(new_input, _inputs[i]);
            }
            return cov_with_inputs;
        }

        const Cov covariance;

        //State
        SVD cov_svd;
        Vector weights;
        std::vector<T> _inputs;
        std::vector<double> _inputs_devs;
    };

}