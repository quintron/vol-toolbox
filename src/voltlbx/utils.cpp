#include "utils.h"

#include <stdexcept>

namespace voltlbx
{
    std::vector<double> linspace(double start, double end, std::size_t size, bool end_point)
    {
        if(size < 2)
            throw std::domain_error("linspace : size should be at least 2");

        if(start >= end)
            throw std::invalid_argument("linspace : start should be strictly lower than end");

        const double step = (end - start) / static_cast<double>(size - 1);
        
        std::vector<double> grid;
        grid.reserve(size);
        for (std::size_t i = 0; i + 1 < size; ++i)
        {
            const double x = start + step * static_cast<double>(i);
            grid.push_back(x);
        }

        if (end_point)
        {
            grid.push_back(end);
        }

        return grid;
    }
}