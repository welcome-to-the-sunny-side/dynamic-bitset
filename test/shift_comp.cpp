#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <random>
#include "dynamic_bitset.hpp"
#include "lazy_dynamic_bitset.hpp"

using dbitset = dybi::dynamic_bitset<uint64_t, 64>;
using ldbitset = dybi::lazy_dynamic_bitset<uint64_t, 64>;

int rng(int n = 1000000007)
{
    static std::mt19937 mt(std::random_device{}());
    return std::uniform_int_distribution<int>(0, n - 1)(mt);
}

const int S = 100;
const int T = 5;

int main()
{
    // int n, q;
    // std::cin >> n >> q;

    constexpr int n = 10000, q = 10000;

    std::vector<bool> v(n);
    std::vector<std::pair<int, int>> op(q);

    for(int i = 0; i < n; i ++)
        v[i] = rng(2);
    
    for(int i = 0; i < q; i ++)
        op[i] = {rng(2), rng(S)};

    int64_t total_a_time = 0, total_b_time = 0, total_c_time = 0;
    int64_t a_tries = 0, b_tries = 0, c_tries = 0;
    
    dbitset a(n);
    ldbitset b(n);
    std::bitset<n> c;

    for(int t = 0; t < 3 * T; t ++)
    {
        if(t % 3 == 0)
        {
            for(int i = 0; i < n; i ++)
                a.set(i, v[i]);

            auto start = std::chrono::high_resolution_clock::now();
            for(int i = 0; i < q; i ++)
            {
                if(op[i].first)
                    a <<= op[i].second;
                else
                    a >>= op[i].second;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_a_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            a_tries ++;
        }
        else if(t % 3 == 1)
        {
            for(int i = 0; i < n; i ++)
                b.set(i, v[i]);

            auto start = std::chrono::high_resolution_clock::now();
            for(int i = 0; i < q; i ++)
            {
                if(op[i].first)
                    b <<= op[i].second;
                else
                    b >>= op[i].second;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_b_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            b_tries ++;
        }
        else
        {
            for(int i = 0; i < n; i ++)
                c[i] = v[i];

            auto start = std::chrono::high_resolution_clock::now();
            for(int i = 0; i < q; i ++)
            {
                if(op[i].first)
                    c <<= op[i].second;
                else
                    c >>= op[i].second;
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_c_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            c_tries ++;
        }
    }
    
    std::cout << "average shift time for regular: " << (long double)total_a_time / (long double)(a_tries * q) << " ns" << std::endl;
    std::cout << "average shift time for lazy: " << (long double)total_b_time / (long double)(b_tries * q) << " ns" << std::endl;
    std::cout << "average shift time for stl: " << (long double)total_c_time / (long double)(c_tries * q) << " ns" << std::endl;
}
