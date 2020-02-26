#pragma once

#include <vector>

namespace voltlbx
{
    /// <summary>
    /// Locate index i such that x is in right open interval [ xs[i], xs[i+1] [ ,  
    /// return -1 if x < xs[0] 
    /// </summary>
    template<typename T>
    int locate_left_index(const std::vector<T>& xs, T x)
    {
        assert(!xs.empty());

        if (x < xs.front())
            return -1;

        if (x >= xs.back())
            return xs.size() - 1;

        int lower = 0;
        int upper = xs.size() - 1;
        while (upper - lower > 1)
        {
            int mid = (upper + lower) >> 1;
            if (x < xs[mid])
            {
                upper = mid;
            }
            else
            {
                lower = mid;
            }
        }
        return lower;
    }
}