#include <iostream>
#include <random>
#include <chrono>
#include "dynamic_bitset.hpp"
#include "lazy_dynamic_bitset.hpp"

using dbitset = dybi::dynamic_bitset;
using ldbitset = dybi::lazy_dynamic_bitset;

int rng(int n = 1000000007)
{
    static std::mt19937 mt(std::random_device{}());
    return std::uniform_int_distribution<int>(0, n - 1)(mt);
}

const int S = 1000;
const int T = 10;

int main()
{
    // int n, q;
    // std::cin >> n >> q;

    constexpr int n = 100000, q = 100000;

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

    int f = 0;

    for(int t = 0; t < T; t ++)
    {
        int idx = rng(n);
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
            f += a[idx];
        }
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
            f += b[idx];
        }
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
            f += c[idx];
        }
    }

    assert(f != -1);
    
    std::cout << "average shift update time for regular dbitset: " << (long double)total_a_time / (long double)(a_tries * q) << " ns" << std::endl;
    std::cout << "average shift update time for lazy dbitset: " << (long double)total_b_time / (long double)(b_tries * q) << " ns" << std::endl;
    std::cout << "average shift update time for stl bitset: " << (long double)total_c_time / (long double)(c_tries * q) << " ns" << std::endl;
}
