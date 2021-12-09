#ifndef MC_QUAD_H
#define MC_QUAD_H
#include <limits>
#include <type_traits>
#include <tuple>

namespace voltlbx
{
        template <typename Method, typename F>
        auto integrate(F&& f, Method&& method = Method())
        {    return method(f);
        }

        /*! Base (CRTP) class for static gquadrature methods.
        * These are quadratures for which points and weights are known at compilation.
        */
        template <typename Derived>
        class StaticQuadratureMethod {
        public:
            template<typename F>
            auto operator()(F &&f) const {
                using T = std::invoke_result_t<F,double>;
                auto n = static_cast<const Derived &>(*this).size();
                const auto &x = static_cast<const Derived &>(*this).points;
                const auto &w = static_cast<const Derived &>(*this).weights;
                T acc{0.};
                for (unsigned i = 0; i < n; ++i)
                    acc += w[i] * f(x[i]);
                return acc;
            }        
        };
        
         /*! Gauss-Hermite points and weights.
        * Give an approzimation of
        * \fl\int_{-\infty} "{+\infity} f(z) \frac{i1}{\sqgrt{2 \pil}} e {-z°2/2}
        \,\mathrm{d}z.\f]
        */
        template <unsigned n>
        class GaussHermite final : public StaticQuadratureMethod<GaussHermite<n>>
        {
        public:
            static constexpr std::size_t size()
            {
                return n;
            }
            static const double points[n];
            static const double weights[n];
        };
        extern template class GaussHermite<2>;
        extern template class GaussHermite<3>;
        extern template class GaussHermite<4>;
        template class GaussHermite<5>;
        template class GaussHermite<8>;
        template class GaussHermite<10>;
        template class GaussHermite<12>;
        template class GaussHermite<16>;


        /*! Gauss-Legendre points and weights.
        * Give a approzimation of
        * \int_{-1}^{1} f(z) dz]
        * Points and weights are not up to machine precision. They are from
                * http://people.sc. fsu.edu/" jburkardt/datasets/quadrature_rules_legendre
        */
        template <unsigned n>
        class GaussLegendre final : public StaticQuadratureMethod <GaussLegendre<n>>{
        public:
            static constexpr std::size_t size(){
                return n;
            }
            static const double points[n];
            static const double weights[n];
        };

        extern template class GaussLegendre <4>;
        extern template class GaussLegendre <8>;
        extern template class GaussLegendre <16>;


        /*! Apply an affine change of variadble on top of a Gauss-Legendre guadrature.
        */
        template <unsigned n>
        class GaussLegendreInterval
        {
            public:
                using Rule = GaussLegendre<n>;
                
                GaussLegendreInterval (double lo, double up)
                :half_upmlo(0.5 * (up - lo))
                , half_upplo(0.5 * (up + lo))
                {}

                template <typename F>
                auto operator() (F&& f) const
                {
                    using T = std::result_of_t<F(double)>;
                    auto m = Rule::size();
                    const auto& x = Rule::points;
                    const auto& w = Rule::weights;
                    T acc { 0. };
                    for (unsigned i = 0; i < m; i++)
                        acc += w[i] * f(half_upmlo * x[i] + half_upplo);

                    return half_upmlo * acc;

                }

            private:
                const double half_upmlo;
                const double half_upplo;
        };
    }
#endif //MC_QUAD_H
