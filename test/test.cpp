#include <iostream>
#include <random>
#include <chrono>
#include "dynamic_bitset.hpp"
#include "lazy_dynamic_bitset.hpp"

using namespace dybi;

int rng(int n = 1000000007)
{
    static std::mt19937 mt(std::random_device{}());
    return std::uniform_int_distribution<int>(0, n - 1)(mt);
}

int main()
{
    dynamic_bitset a(100000, true);

    std::vector<dynamic_bitset> store (1000, dynamic_bitset(100000, true));

    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < 10000; i ++)
    {
        int type = rng(3);
        int idx = rng(store.size());

        if(type == 0)
            a |= store[idx];
        else if(type == 1)
            a &= store[idx];
        else
            a ^= store[idx];
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << "Time: " << time << " ns" << std::endl;
}