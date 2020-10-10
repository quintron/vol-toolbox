#include <gtest/gtest.h>
#include <voltlbx/linear_filter.h>
#include <voltlbx/utils.h>

using namespace voltlbx;

TEST(LinearFilter, RankOne)
{
    class RankOneCov : public Covariance<double>
    {
    public:

        RankOneCov(double dev)
            :_dev(dev)
        {
        }

        double operator()(double, double) const override
        {
            return _dev * _dev;
        }

        double dev(double) const override
        {
            return _dev;
        }

        double correl(double, double) const override
        {
            return 1.0;
        }

    private:
        const double _dev;
    };

    const double dev = 0.01;
    auto filter = LinearFilter<double>(std::make_shared<RankOneCov>(dev));
    
    const std::vector<double> xs(50, NAN);
    
    const double y = 3.1415;
    const std::vector<double> ys(50, y);

    const double eps = 0.001;
    const std::vector<double> err_devs(50, eps);

    filter.fit(xs, ys, err_devs);

    const double ref_pred = y / (1.0 + eps * eps / (dev * dev * ys.size()));
    for (auto x : linspace(-10.0, 10.0, 50))
    {
        auto p = filter.predict(x);
        ASSERT_NEAR(p, ref_pred, 1.0e-15);
    }
}
