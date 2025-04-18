#include <vector>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <bitset>

namespace dybi
{
    constexpr int AVCT = 16;

    template<typename T, const int B>
    class dynamic_bitset
    {
        public:

        static_assert(sizeof(T) * 8 == B, "check block width");
        static_assert(std::is_same<T, uint64_t>::value, "modify popcnt(), ctz(), clz()");

        static inline constexpr bool on(int i, T x) noexcept
        {
            return ((T(1) << i) & x) != 0;
        }
        static inline constexpr T prefix(int i) noexcept
        {
            return (i >= B) ? ~T(0) : ((T(1) << i) - T(1));
        }
        static inline constexpr T suffix(int i) noexcept
        {
            return ~prefix(B - i);
        }
        static inline constexpr T range(int l, int r) noexcept
        {
            return prefix(r) ^ prefix(l - 1);
        }
        static constexpr int popcnt(T x) noexcept
        {
            // return _mm_popcnt_u64(x);
            return __builtin_popcountll(x);
        }
        static constexpr int clz(T x) noexcept
        {
            return __builtin_clzll(x);
        }
        static constexpr int ctz(T x) noexcept
        {
            return __builtin_ctzll(x);
        }
        static inline constexpr int block_id(int i) noexcept
        {
            return i / B;
        }

        int n, m;
        std::vector<T> b;

        // helper functions

        // returns a submask of a single block 
        inline T submask(int l, int r) const noexcept
        {
            int bx = block_id(l);
            assert(bx == block_id(r));
            return (b[bx] & range(l - bx * B + 1, r - bx * B + 1)); 
        }
        // trim the overhang of the last block 
        inline void trim() noexcept
        {
            b.back() &= prefix(n % B == 0 ? B : n % B);
        }

        dynamic_bitset (int n) : dynamic_bitset(n, false) {};
        dynamic_bitset(int n, bool init) : n(n), m((n + B - 1)/B), b(m, init ? ~T(0) : T(0)) 
        {
            trim();
        };

        // set the i-th bit to val
        inline void set(int i, bool val) noexcept
        {
            assert(0 <= i and i < n);
            if(val)
                b[i/B] |= (T(1) << (i % B));
            else
                b[i/B] &= ~(T(1) << (i % B));
        }

        // get the value of the i-th bit
        inline bool get(int i) const noexcept
        {
            assert(0 <= i and i < n);
            return (b[i/B] & (T(1) << (i % B))) != 0;
        }

        // reset the bitset
        void reset() noexcept
        {
            std::fill(b.begin(), b.end(), T(0));
        }

        // bitwise operations
        
        // note : if the other bitset is smaller - it is padded with zeros
        // if the other bitset is larger - the overhanging suffix is ignored

        void operator &= (const dynamic_bitset &other)
        {
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] &= other.b[i];
            if(m > other.m)
                std::fill(b.begin() + other.m, b.begin() + m, T(0));
            // trim(); there's no need to trim here, as no extra bits are switched on in our own overhang
        }

        void operator |= (const dynamic_bitset &other)
        {
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] |= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
    
        void operator ^= (const dynamic_bitset &other)
        {
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] ^= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }

        // shift left by x bits
        void operator <<= (int x)
        {
            if(x == 0)
                return;
    
            if(x >= n)
            {
                reset();
                return;
            }

            const int s = x/B, d = x % B, r = B - d;
    
            if(d > 0)
            {
                // #pragma GCC unroll AVCT
                #pragma ivdep
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = (b[i] << d) | (b[i - 1] >> r);
                b[s] = b[0] << d;
            }
            else
            {
                // #pragma GCC unroll AVCT
                #pragma ivdep
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = b[i];
                b[s] = b[0];
            }

            std::fill(b.begin(), b.begin() + s, T(0));
    
            trim();
        }

        // shift right by x bits
        void operator >>= (int x)
        {
            if(x == 0)
                return;
        
            if(x >= n)
            {
                reset();
                return;
            }
    
            const int s = x/B, d = x % B, l = B - d;
    
            if(d > 0)
            {
                // #pragma GCC unroll AVCT
                #pragma ivdep
                for(int i = s; i < m - 1; i ++)
                    b[i - s] = (b[i] >> d) | (b[i + 1] << l); 
                b[m - 1 - s] = b[m - 1] >> d;
            }
            else
            {
                // #pragma GCC unroll AVCT
                #pragma ivdep
                for(int i = s; i < m; i ++)
                    b[i - s] = b[i];
            }
    
            std::fill(b.begin() + m - s, b.end(), T(0));        
    
            // trim();
        }

        // equality requires equal size and contents
        bool operator == (const dynamic_bitset &other)
        {
            return ((n == other.n) and b == other.b); 
        }

        bool operator != (const dynamic_bitset &other)
        {
            return !(*this == other);
        }


        // more bitwise

        dynamic_bitset operator & (const dynamic_bitset &other)
        {
            dynamic_bitset result(*this);
            result &= other;
            return result;
        }
    
        dynamic_bitset operator | (const dynamic_bitset &other)
        {
            dynamic_bitset result(*this);
            result |= other;
            return result;
        }
    
        dynamic_bitset operator ^ (const dynamic_bitset &other)
        {
            dynamic_bitset result(*this);
            result ^= other;
            return result;
        }
    
        dynamic_bitset operator >> (int x)
        {
            dynamic_bitset result(*this);
            result >>= x;
            return result;
        }
    
        dynamic_bitset operator << (int x)
        {
            dynamic_bitset result(*this);
            result <<= x;
            return result;
        }
    
        dynamic_bitset operator ~()
        {
            dynamic_bitset result(*this);
            for(auto &v : result)
                v = ~v;
            result.trim();
            return result;
        }


        //custom operations

        // returns the number of set bits
        int count() const noexcept
        {
            return std::accumulate(b.begin(), b.end(), 0, [](int sum, T value) { return sum + popcnt(value); });
        }
        
        // returns the index of the first set bit (-1 if none)
        int find_first()
        {
            int pos = -1;

            for(int bi = 0; bi < m; bi ++)
            {
                if(b[bi] == T(0))
                    continue;
                
                pos = ctz(b[bi]) + bi * B;
                break;
            }

            return pos;
        }

        // returns the index of the last set bit (-1 if none)
        int find_last()
        {
            int pos = -1;

            for(int bi = m - 1; bi >= 0; bi --)
            {
                if(b[bi] == T(0))
                    continue;
                
                pos = B - clz(b[bi]) - 1 + bi * B;
                break;
            }

            return pos;
        }

        // perform an arbitrary operation on the range [l, r]
        // usage: pass two lambdas, `block_brute(l, r)` and `block_quick(block_id)`
        // `block_brute(l, r)` shall be called for each block that lies partially in the range [l, r]
        // `block_quick(block_id)` shall be called for each block that is entirely within [l, r]
        template<typename F1, typename F2>
        void range_process(int l, int r, F1 block_brute, F2 block_quick)
        {
            assert(0 <= l and l <= r and r < n);
    
            int bl = block_id(l), br = block_id(r);
    
            if(bl == br)
                block_brute(l, r);
            else
            {
                block_brute(l, (bl + 1) * B - 1);
                for(int bi = bl + 1; bi < br; bi ++)
                    block_quick(bi);
                block_brute(br * B, r);
            }
        }
    
        // some helpful `range_process` clients

        // set the range [l, r] to val
        void range_set(int l, int r, bool val)
        {
            auto block_brute = [&](int l, int r) -> void
            {
                int bi = block_id(l);
                T mask = range(l - bi * B + 1, r - bi * B + 1);
                if(val)
                    b[bi] |= mask;
                else
                    b[bi] &= ~mask;
            };
            auto block_quick = [&](int bi) -> void
            {
                b[bi] = (val ? ~T(0) : T(0));
            };
            range_process(l, r, block_brute, block_quick);
        }

        // count the number of set bits in the range [l, r]
        int count(int l, int r)
        {
            int cnt = 0;
            auto block_brute = [&](int l, int r) -> void
            {
                cnt += popcnt(submask(l, r));
            };
            auto block_quick = [&](int bi) -> void
            {
                cnt += popcnt(b[bi]);
            };
            range_process(l, r, block_brute, block_quick);
            return cnt;
        }

        // returns the index of the first set bit in the range [l, r] (-1 if none)
        int find_first (int l, int r)
        {
            int pos = -1;
            auto block_brute = [&](int l, int r) -> void
            {
                for(int i = l; i <= r and pos == -1; i ++)
                    if(get(i))
                        pos = i;    
            };
            auto block_quick = [&](int bi) -> void
            {
                if(b[bi] == T(0) or pos != -1)
                    return;
    
                pos = ctz(b[bi]) + bi * B;
            };
    
            range_process(l, r, block_brute, block_quick);
            return pos;
        }

        // returns the index of the last set bit in the range [l, r] (-1 if none)
        int find_last(int l, int r)
        {
            int pos = -1;
            auto block_brute = [&](int l, int r) -> void
            {
                for(int i = l; i <= r; i ++)
                    if(get(i))
                        pos = i;    
            };
            auto block_quick = [&](int bi) -> void
            {
                if(b[bi] == T(0))
                    return;
    
                pos = B - clz(b[bi]) - 1 + bi * B;
            };
    
            range_process(l, r, block_brute, block_quick);
            return pos;
        }

        friend std::ostream &operator<<(std::ostream &os, const dynamic_bitset &bitset)
        {
            for (int i = bitset.m - 1; i >= 0; --i)
                os << std::bitset<B>(bitset.b[i]);
            os << '\n';
            return os;
        }
    };
}
