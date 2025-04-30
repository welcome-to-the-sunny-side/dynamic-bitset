#include <vector>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <bitset>

namespace dybi
{
    constexpr int L_AVCT = 16;
    class lazy_dynamic_bitset
    {
        public:
        static constexpr int B = 64;

        class bit_reference
        {
            lazy_dynamic_bitset &bitset_ref;
            int pos;

        public:
            bit_reference(lazy_dynamic_bitset &ref, int p) : bitset_ref(ref), pos(p) {}

            // Assignment operator for writing (e.g., bitset[i] = true)
            bit_reference &operator=(bool val)
            {
                bitset_ref.resolve_shift(); // Resolve before modifying
                bitset_ref.set(pos, val);
                return *this;
            }

            // Assignment operator for copying from another reference (e.g., bitset[i] = bitset[j])
            bit_reference &operator=(const bit_reference &other)
            {
                // Need to resolve both if the other one could be lazy too,
                // but 'other' is const, so we cast away constness only to call resolve_shift.
                // It's generally safe as resolve_shift only modifies mutable members.
                const_cast<lazy_dynamic_bitset &>(other.bitset_ref).resolve_shift();
                bitset_ref.resolve_shift();
                bitset_ref.set(pos, static_cast<bool>(other)); // Reading other calls its operator bool()
                return *this;
            }

            // Conversion operator for reading (e.g., if (bitset[i]))
            operator bool() const
            {
                const_cast<lazy_dynamic_bitset &>(bitset_ref).resolve_shift(); // Resolve before reading
                return bitset_ref.get(pos);
            }

            // Flip the bit (e.g., bitset[i].flip())
            bit_reference &flip()
            {
                bitset_ref.resolve_shift(); // Resolve before modifying
                bitset_ref.set(pos, !bitset_ref.get(pos));
                return *this;
            }

            bool operator~() const
            {
                const_cast<lazy_dynamic_bitset &>(bitset_ref).resolve_shift(); // Resolve before reading
                return !bitset_ref.get(pos);
            }
        };

        static_assert(sizeof(uint64_t) * 8 == B, "check block width");
        static_assert(std::is_same<uint64_t, uint64_t>::value, "modify popcnt(), ctz(), clz()");

        static inline constexpr bool on(int i, uint64_t x) noexcept
        {
            return ((uint64_t(1) << i) & x) != 0;
        }
        static inline constexpr uint64_t prefix(int i) noexcept
        {
            return (i >= B) ? ~uint64_t(0) : ((uint64_t(1) << i) - uint64_t(1));
        }
        static inline constexpr uint64_t suffix(int i) noexcept
        {
            return ~prefix(B - i);
        }
        static inline constexpr uint64_t range(int l, int r) noexcept
        {
            return prefix(r) ^ prefix(l - 1);
        }
        static constexpr int popcnt(uint64_t x) noexcept
        {
            // return _mm_popcnt_u64(x);
            return __builtin_popcountll(x);
        }
        static constexpr int clz(uint64_t x) noexcept
        {
            return __builtin_clzll(x);
        }
        static constexpr int ctz(uint64_t x) noexcept
        {
            return __builtin_ctzll(x);
        }
        static inline constexpr int block_id(int i) noexcept
        {
            return i / B;
        }

        int n, m;
        mutable std::vector<uint64_t> b;
        mutable int pending_shift;           //positive => left shift, negative => right shift

        // helper functions

        // returns a submask of a single block 
        inline uint64_t submask(int l, int r) const noexcept
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
        lazy_dynamic_bitset(int n, bool init) : n(n), m((n + B - 1)/B), b(m, init ? ~uint64_t(0) : uint64_t(0)), pending_shift(0) 
        {
            trim();
        };

        // set the i-th bit to val
        inline void set(int i, bool val) noexcept
        {
            assert(0 <= i and i < n);
            resolve_shift();
            if(val)
                b[i/B] |= (uint64_t(1) << (i % B));
            else
                b[i/B] &= ~(uint64_t(1) << (i % B));
        }

        // get the value of the i-th bit
        inline bool get(int i) const noexcept
        {
            assert(0 <= i and i < n);
            resolve_shift();
            return (b[i/B] & (uint64_t(1) << (i % B))) != 0;
        }

        // reset the bitset
        void reset() noexcept
        {
            std::fill(b.begin(), b.end(), uint64_t(0));
        }

        // bitwise operations
        
        // note : if the other bitset is smaller - it is padded with zeros
        // if the other bitset is larger - the overhanging suffix is ignored

        void operator &= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            #pragma ivdep
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] &= other.b[i];
            if(m > other.m)
                std::fill(b.begin() + other.m, b.begin() + m, uint64_t(0));
            // trim(); there's no need to trim here, as no extra bits are switched on in our own overhang
        }

        void operator |= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            #pragma ivdep
            for(int i = 0; i < std::min(m, other.m); i ++)
                b[i] |= other.b[i];
            trim(); // this might result in some overhanging bits being switched on
        }
    
        void operator ^= (const lazy_dynamic_bitset &other)
        {
            resolve_shift();
            other.resolve_shift();

            #pragma ivdep
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
                #pragma ivdep
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = (b[i] << d) | (b[i - 1] >> r);
                b[s] = b[0] << d;
            }
            else
            {
                #pragma ivdep
                for(int i = m - 1 - s; i > 0; i --)
                    b[i + s] = b[i];
                b[s] = b[0];
            }
    
            std::fill(b.begin(), b.begin() + s, uint64_t(0));
    
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
                #pragma ivdep
                for(int i = s; i < m - 1; i ++)
                    b[i - s] = (b[i] >> d) | (b[i + 1] << l); 
                b[m - 1 - s] = b[m - 1] >> d;
            }
            else
            {
                #pragma ivdep
                for(int i = s; i < m; i ++)
                    b[i - s] = b[i];
            }
    
            std::fill(b.begin() + m - s, b.end(), uint64_t(0));        
    
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
            for(auto &v : result.b)
                v = ~v;
            result.trim();
            return result;
        }

        //custom operations

        // returns the number of set bits
        int count() const noexcept
        {
            resolve_shift();
            return std::accumulate(b.begin(), b.end(), 0, [](int sum, uint64_t value) { return sum + popcnt(value); });
        }
        
        // returns the index of the first set bit (-1 if none)
        int find_first()
        {
            resolve_shift();
            int pos = -1;

            for(int bi = 0; bi < m; bi ++)
            {
                if(b[bi] == uint64_t(0))
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
                if(b[bi] == uint64_t(0))
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
                uint64_t mask = range(l - bi * B + 1, r - bi * B + 1);
                if(val)
                    b[bi] |= mask;
                else
                    b[bi] &= ~mask;
            };
            auto block_quick = [&](int bi) -> void
            {
                b[bi] = (val ? ~uint64_t(0) : uint64_t(0));
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
                if(b[bi] == uint64_t(0) or pos != -1)
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
                if(b[bi] == uint64_t(0))
                    return;
    
                pos = B - clz(b[bi]) - 1 + bi * B;
            };
    
            range_process(l, r, block_brute, block_quick);
            return pos;
        }

        // Non-const version: returns a modifiable reference proxy
        bit_reference operator[](int i)
        {
            assert(0 <= i and i < n);
            // Don't resolve shift here, the reference object will do it on access/modification
            return bit_reference(*this, i);
        }

        // Const version: returns the actual bool value (read-only)
        bool operator[](int i) const
        {
            assert(0 <= i and i < n);
            resolve_shift(); // Resolve shift before getting
            return get(i);   // Directly call get for const access
        }

        friend std::ostream &operator<<(std::ostream &os, const lazy_dynamic_bitset &bitset)
        {
            // No need to resolve shift explicitly here, operator[] const will do it if needed
            // or the internal get used by the reference proxy will handle it.
            // However, printing the whole bitset requires resolving.
            bitset.resolve_shift();
            for (int i = bitset.m - 1; i >= 0; --i)
                os << std::bitset<B>(bitset.b[i]);
            os << '\n';
            return os;
        }
    };
}
