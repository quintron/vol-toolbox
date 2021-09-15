#include "optim.h"

#include <cminpack.h>

namespace voltlbx
{
    int MinPackWrapper::solve(std::size_t n, double * x, std::size_t m, double * f, double tol)
    {
        auto v_iwa = std::vector<int>(n);
        auto v_wa = std::vector<double>(n*m + 5 * n + m);
        auto nn = static_cast<int>(n);
        auto mm = static_cast<int>(m);
        auto pp = static_cast<int>(v_wa.size());
        auto iwa = v_iwa.data();
        auto wa = v_wa.data();

        int res = lmdif1(&wrapper, this, mm, nn, x, f, tol, iwa, wa, pp);
        switch (res)
        {
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
            return res;

        default:
        case -1:
            if (error)
                std::rethrow_exception(error);
            else
                throw std::runtime_error("execution error in levmarq");
        case 0:
            throw std::logic_error("improper input parameters in levmarq");
        }
    }

    int MinPackWrapper::wrapper(void *p, int m, int n, const double * x, double * f, int i)
    {
        auto self = static_cast<MinPackWrapper *>(p);

        try 
        {
            auto nn = static_cast<std:: size_t>(n);
            auto mm = static_cast<std::size_t>(m);
            self->func(nn, x, mm, f);
            return i;
        }
        catch (...)
        {
            self->error = std::current_exception();
            return -1;
        }        
    }

}