// #pragma GCC optimize("O3,unroll-loops")
// if available, enable AVX2
// #pragma GCC target("avx2")

#include <iostream>
#include <random>
#include "dynamic_bitset.hpp"

using dbitset = dybi::dynamic_bitset<uint64_t, 64>;

int rng(int n = 1000000007)
{
    static std::mt19937 mt(std::random_device{}());
    return std::uniform_int_distribution<int>(0, n - 1)(mt);
}

int main()
{
    int n, q;
    std::cin >> n >> q;

    dbitset a(n);

    for(int i = 0; i < n; i ++)
    {
        int x;
        std::cin >> x;
        a.set(i, x);
    }

    for(int ts = 0; ts < q; ts ++)
    {
        int t;
        std::cin >> t;
        if(t == 1)
        {
            int i, x;
            std::cin >> i >> x;
            a.set(i, x);
        }
        else if(t == 2)
        {
            int i;
            std::cin >> i;
            std::cout << a.get(i) << std::endl;
        }
        else if(t == 3)
        {
            int l, r;
            std::cin >> l >> r;
            std::cout << a.count(l, r) << std::endl;
        }
        else if(t == 4)
        {
            int x;
            std::cin >> x;
            a >>= x;
        }
        else if(t == 5)
        {
            int x;
            std::cin >> x;
            a <<= x;
        }
    }
}