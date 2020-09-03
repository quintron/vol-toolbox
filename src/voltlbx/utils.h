#pragma once

#include <algorithm>
#include <vector>
#include <map>
#include <type_traits>

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
            return static_cast<int>(xs.size()) - 1;

        int lower = 0;
        int upper = static_cast<int>(xs.size()) - 1;
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


    /*! Return the longest sequence of integer a[i] s.t. : 
        1. 0 <= a[i] < N
        2. pred(a[i], a[i+1]) is true
    */
    template<typename BinaryPredicate>
    std::vector<std::size_t> longest_chain(std::size_t N, BinaryPredicate&& pred)
    {       
        if (N == 0)
            return {};
        
        // ls[i] will contain the size of longest chain with endpoint a[i]
        std::vector<std::size_t> ls(N, 1);
        for (std::size_t i = 1; i < N; ++i)
        {
            for (std::size_t j = 0; j < i; ++j)
            {
                if (ls[i] < ls[j] + 1
                    && pred(j, i))
                {
                    ls[i] = ls[j] + 1;
                }
            }                                   
        }

        const auto max_it = std::max_element(std::begin(ls), std::end(ls));
        const std::size_t max_index = std::distance(std::begin(ls), max_it);

        //Backtracking            
        std::size_t current_size = *max_it;
        std::size_t current_index = max_index;
        std::vector<std::size_t> chain = { max_index };
        for (int i = static_cast<int>(max_index) - 1; i >= 0; --i)
        {
            if (ls[i] + 1 == current_size
                && pred(i, current_index))
            {
                --current_size;
                current_index = static_cast<std::size_t>(i);
                chain.push_back(current_index);                
            }
        }
        std::sort(std::begin(chain), std::end(chain));
        return chain;
    }


    template<typename T>
    std::vector<std::size_t> longest_increasing_subsequence(const std::vector<T>& a)
    {
        return longest_chain(std::size(a), [&](std::size_t i, std::size_t j) {return a[i] < a[j]; });
    }


    std::vector<double> linspace(double start, double end, std::size_t size, bool end_point = true);

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template<typename Tin, typename F>
    auto map(const std::vector<Tin>& in, F&& func)
    {
        using Tout = std::invoke_result<F, Tin>::type;
        std::vector<Tout> res(in.size());
        std::transform(in.begin(), in.end(), res.begin(), func);
        return res;
    }

}