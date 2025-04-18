#include <bits/stdc++.h>


namespace dybi
{
    template<typename T, const int B>
    class lazy_dynamic_bitset
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
        mutable std::vector<T> b;
        mutable int pending_shift;           //positive => left shift, negative => right shift

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

        lazy_dynamic_bitset (int n) : lazy_dynamic_bitset(n, false) {};
        lazy_dynamic_bitset(int n, bool init) : n(n), m((n + B - 1)/B), b(m, init ? ~T(0) : T(0)), pending_shift(0) 
        {
            trim();
        };

        // set the i-th bit to val
        inline void set(int i, bool val) noexcept
        {
            assert(0 <= i and i < n);
            resolve_shift();
            if(val)
                b[i/B] |= (T(1) << (i % B));
            else
                b[i/B] &= ~(T(1) << (i % B));
        }

        // get the value of the i-th bit
        inline bool get(int i) const noexcept
        {
            assert(0 <= i and i < n);
            resolve_shift();
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

        void operator &= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] &= other.b[i];
            if(m > other.m)
                std::fill(b.begin() + other.m, b.begin() + m, T(0));
            // trim(); there's no need to trim here, as no extra bits are switched on in our own overhang
        }

        void operator |= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] |= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
    
        void operator ^= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] ^= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }


        // shift operators

        // immediate shifters
        void left_shift(int x)
        {
            if(x >= n)
            {
                reset();
                return;
            }

            const int s = x/B, d = x % B, r = B - d;
    
            if(d > 0)
            {
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = (b[i] << d) | (b[i - 1] >> r);
                b[s] = b[0] << d;
            }
            else
            {
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = b[i];
                b[s] = b[0];
            }
    
            std::fill(b.begin(), b.begin() + s, T(0));
    
            trim();
        }

        // immediate right shift
        void right_shift(int x)
        {
            if(x >= n)
            {
                reset();
                return;
            }
    
            const int s = x/B, d = x % B, l = B - d;
    
            if(d > 0)
            {
                for(int i = s; i < m - 1; i ++)
                    b[i - s] = (b[i] >> d) | (b[i + 1] << l); 
                b[m - 1 - s] = b[m - 1] >> d;
            }
            else
                for(int i = s; i < m; i ++)
                    b[i - s] = b[i];
    
            std::fill(b.begin() + m - s, b.end(), T(0));        
    
            // trim();
        }

        inline void resolve_shift() const noexcept
        {
            if(pending_shift == 0)
                return;
            auto self = const_cast<lazy_dynamic_bitset*>(this);   // safe: we only
            if(pending_shift > 0)
                self->left_shift(pending_shift);
            else
                self->right_shift(-pending_shift);
            pending_shift = 0;
        }

        // lazily shift left by x bits (if LAZY = true)
        void operator <<= (int x)
        {
            if(pending_shift < 0)
                resolve_shift();
            pending_shift += x;
            return;
        }

        // lazily shift right by x bits (if LAZY = true)
        void operator >>= (int x)
        {
            if(pending_shift > 0)
                resolve_shift();
            pending_shift -= x;
            return;
        }

        // equality requires equal size and contents
        bool operator == (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();
            return ((n == other.n) and b == other.b); 
        }

        bool operator != (const lazy_dynamic_bitset &other)
        {
            return !(*this == other);
        }

        // more bitwise

        lazy_dynamic_bitset operator & (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            lazy_dynamic_bitset result(*this);
            result &= other;
            return result;
        }
    
        lazy_dynamic_bitset operator | (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            lazy_dynamic_bitset result(*this);
            result |= other;
            return result;
        }
    
        lazy_dynamic_bitset operator ^ (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            lazy_dynamic_bitset result(*this);
            result ^= other;
            return result;
        }
    
        lazy_dynamic_bitset operator >> (int x)
        {
            lazy_dynamic_bitset result(*this);
            result.right_shift(x);          //avoid lazy evaluation
            return result;
        }
    
        lazy_dynamic_bitset operator << (int x)
        {
            lazy_dynamic_bitset result(*this);
            result.left_shift(x);          //avoid lazy evaluation
            return result;
        }
    
        lazy_dynamic_bitset operator ~()
        {
            resolve_shift();
            lazy_dynamic_bitset result(*this);
            for(auto &v : result)
                v = ~v;
            result.trim();
            return result;
        }

        //custom operations

        // returns the number of set bits
        int count() const noexcept
        {
            resolve_shift();
            return std::accumulate(b.begin(), b.end(), 0, [](int sum, T value) { return sum + popcnt(value); });
        }
        
        // returns the index of the first set bit (-1 if none)
        int find_first()
        {
            resolve_shift();
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
            resolve_shift();
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
            resolve_shift();

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

        friend std::ostream &operator<<(std::ostream &os, const lazy_dynamic_bitset &bitset)
        {
            bitset.resolve_shift();
            for (int i = bitset.m - 1; i >= 0; --i)
                os << std::bitset<B>(bitset.b[i]);
            os << '\n';
            return os;
        }
    };
}
namespace dybi
{
    #define AVCT 16
    #define SIMD_EXPLICIT 0

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

        #if defined(SIMD_EXPLICIT) && defined(__AVX2__)
        void operator &= (const dynamic_bitset &other)
        {
            constexpr int VEC_SIZE = sizeof(__m256i) / sizeof(T); // Number of T elements in a 256-bit vector
            const int limit = std::min(m, other.m);
            const int vec_end = limit - (limit % VEC_SIZE);

            // Process full 256-bit chunks using AVX2
            for (int i = 0; i < vec_end; i += VEC_SIZE) 
            {
                // Load 256 bits (4 x uint64_t) from both bitsets
                __m256i vec_a = _mm256_loadu_si256((__m256i*)&b[i]);
                __m256i vec_b = _mm256_loadu_si256((__m256i*)&other.b[i]);

                // Perform AND
                __m256i result_vec = _mm256_and_si256(vec_a, vec_b);

                // Store the result back
                _mm256_storeu_si256((__m256i*)&b[i], result_vec);
            }

            // Process any remaining elements scalar way
            for (int i = vec_end; i < limit; ++i) 
                b[i] &= other.b[i];

            trim(); // this might result in some overhanging bits being switched on
        }
        #else 
        void operator &= (const dynamic_bitset &other)
        {
            // #pragma GCC unroll AVCT
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] &= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
        #endif

        #if defined(SIMD_EXPLICIT) && defined(__AVX2__)
        void operator |= (const dynamic_bitset &other)
        {
            constexpr int VEC_SIZE = sizeof(__m256i) / sizeof(T); // Number of T elements in a 256-bit vector
            const int limit = std::min(m, other.m);
            const int vec_end = limit - (limit % VEC_SIZE);

            // Process full 256-bit chunks using AVX2
            for (int i = 0; i < vec_end; i += VEC_SIZE) 
            {
                // Load 256 bits (4 x uint64_t) from both bitsets
                __m256i vec_a = _mm256_loadu_si256((__m256i*)&b[i]);
                __m256i vec_b = _mm256_loadu_si256((__m256i*)&other.b[i]);

                // Perform OR
                __m256i result_vec = _mm256_or_si256(vec_a, vec_b);

                // Store the result back
                _mm256_storeu_si256((__m256i*)&b[i], result_vec);
            }

            // Process any remaining elements scalar way
            for (int i = vec_end; i < limit; ++i) 
                b[i] |= other.b[i];

            trim(); // this might result in some overhanging bits being switched on
        }
        #else 
        void operator |= (const dynamic_bitset &other)
        {
            // #pragma GCC unroll AVCT
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] |= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
        #endif

        #if defined(SIMD_EXPLICIT) && defined(__AVX2__)
        void operator ^= (const dynamic_bitset &other)
        {
            constexpr int VEC_SIZE = sizeof(__m256i) / sizeof(T); // Number of T elements in a 256-bit vector
            const int limit = std::min(m, other.m);
            const int vec_end = limit - (limit % VEC_SIZE);

            // Process full 256-bit chunks using AVX2
            for (int i = 0; i < vec_end; i += VEC_SIZE) 
            {
                // Load 256 bits (4 x uint64_t) from both bitsets
                __m256i vec_a = _mm256_loadu_si256((__m256i*)&b[i]);
                __m256i vec_b = _mm256_loadu_si256((__m256i*)&other.b[i]);

                // Perform XOR
                __m256i result_vec = _mm256_xor_si256(vec_a, vec_b);

                // Store the result back
                _mm256_storeu_si256((__m256i*)&b[i], result_vec);
            }

            // Process any remaining elements scalar way
            for (int i = vec_end; i < limit; ++i) 
                b[i] ^= other.b[i];

            trim(); // this might result in some overhanging bits being switched on
        }
        #else 
        void operator ^= (const dynamic_bitset &other)
        {
            // #pragma GCC unroll AVCT
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] ^= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
        #endif

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
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = (b[i] << d) | (b[i - 1] >> r);
                b[s] = b[0] << d;
            }
            else
            {
                #pragma GCC unroll AVCT
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
                for(int i = s; i < m - 1; i ++)
                    b[i - s] = (b[i] >> d) | (b[i + 1] << l); 
                b[m - 1 - s] = b[m - 1] >> d;
            }
            else
            {
                // #pragma GCC unroll AVCT
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


using dbitset = dybi::dynamic_bitset<uint64_t, 64>;
using ldbitset = dybi::lazy_dynamic_bitset<uint64_t, 64>;

int rng(int n = 1000000007)
{
    static std::mt19937 mt(std::random_device{}());
    return std::uniform_int_distribution<int>(0, n - 1)(mt);
}

const int P = 100;
const int T = 5;

int main()
{
    // int n, q;
    // std::cin >> n >> q;

    const int n = 10000, q = 10000;
    
    std::vector<dbitset> as (P, dbitset(n));
    std::vector<ldbitset> bs (P, ldbitset(n));
    std::vector<std::bitset<n>> cs (P);

    for(int i = 0; i < P; i ++)
    {
        for(int j = 0; j < n; j ++)
        {
            int x = rng(2);
            as[i].set(j, x);
            bs[i].set(j, x);
            cs[i].set(j, x);
        }
    }

    std::vector<bool> v(n);
    std::vector<std::pair<int, int>> op(q);

    for(int i = 0; i < n; i ++)
        v[i] = rng(2);
    
    for(int i = 0; i < q; i ++)
        op[i] = {rng(3), rng(P)};

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
                a ^= as[op[i].second];
                // if(op[i].first == 0)
                    // a &= (as[op[i].second] | as[op[i].second]);
                // else if(op[i].first == 1)
                    // a |= (as[op[i].second] ^ as[op[i].second]);
                // else if(op[i].first == 2)
                    // a ^= (as[op[i].second] & as[op[i].second]);
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
                b ^= bs[op[i].second];
                // if(op[i].first == 0)
                    // b &= (bs[op[i].second] | bs[op[i].second]);
                // else if(op[i].first == 1)
                    // b |= (bs[op[i].second] ^ bs[op[i].second]);
                // else if(op[i].first == 2)
                    // b ^= (bs[op[i].second] & bs[op[i].second]);
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
                c ^= cs[op[i].second];  
                // if(op[i].first == 0)
                    // c &= (cs[op[i].second] | cs[op[i].second]);
                // else if(op[i].first == 1)
                    // c |= (cs[op[i].second] ^ cs[op[i].second]);
                // else if(op[i].first == 2)
                    // c ^= (cs[op[i].second] & cs[op[i].second]);
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_c_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            c_tries ++;
        }
    }

    long double reg_time = (long double)total_a_time / (long double)(a_tries * q);
    long double lazy_time = (long double)total_b_time / (long double)(b_tries * q);
    long double stl_time = (long double)total_c_time / (long double)(c_tries * q);
    
    // std::cout << "average update time for regular: " << reg_time << " ns" << std::endl;
    // std::cout << "average update time for lazy: " << lazy_time << " ns" << std::endl;
    // std::cout << "average update time for stl: " << (long double)total_c_time / (long double)(c_tries * q) << " ns" << std::endl;
    std::cout << reg_time << "_" << lazy_time << "_" << stl_time << std::endl;
}
