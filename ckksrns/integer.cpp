
#include "integer.h"
#include "err.h"
#include "lifter.h"

Integer modmul(Integer ia, Integer ib, Integer imod)
{
    using u = uint64_t;
    if (sizeof(u) != sizeof(Integer)) nevers("not implemented");

    u a(ia), b(ib), m(imod);
    if (!((a >> 32) | (b >> 32)))
        return Integer(a * b % m);

    using u128 = lifter<u>;
    u128 ua { a }, ub { b }, um { m };
    Integer result = ((ua * ub) % um).lo;
    return result;
}

Integer modadd(Integer ia, Integer ib, Integer imod)
{
    using u = uint64_t;
    if (sizeof(u) != sizeof(Integer)) nevers("not implemented");

    u a(ia), b(ib), m(imod);
    if (!((a >> 63) | (b >> 63)))
        return Integer((a + b) % m);

    using u128 = lifter<u>;
    u128 ua { a }, ub { b }, um { m };
    Integer result = ((ua + ub) % um).lo;
    return result;
}

Integer modsub(Integer ia, Integer ib, Integer imod)
{
    return modadd(ia, imod - ib, imod);
}

template <class T>
T egcdT_rec(const T & a, const T & b, T & s0, T & t0, T s1, T t1, T s2, T t2)
{
    auto c = a % b;
    if (c == T(0)) return b;
    auto q = (a - c) / b; // exact division
    s0 = s2 - s1 * q;
    t0 = t2 - t1 * q;
    return egcdT_rec(b, c, s0, t0, s0, t0, s1, t1);
}

template <class T>
T egcdT(const T & a, const T & b, T & s0, T & t0)
{
    T s1(0), t1(1), s2(1), t2(0);
    return egcdT_rec(a, b, s0, t0, s1, t1, s2, t2);
}

Integer invert(Integer a, Integer mod)
{
    Integer s, t;
    auto gcd = egcdT<Integer>(a, mod, s, t);
    if (gcd * gcd != 1) nevers("no inverse (ull)");
    while (s < 0) s += mod;
    s /= gcd;
    return s % mod;
}

Integer modinv(Integer a, Integer mod)
{
    return invert(a, mod);
}

bool coprime(Integer a, Integer b)
{
    Integer s, t;
    auto gcd = egcdT<Integer>(a, b, s, t);
    return (gcd * gcd == 1);
}

std::pair<Integer, int> modaddOver(Integer ia, Integer ib, Integer imod)
{
    using u = uint64_t;
    if (sizeof(u) != sizeof(Integer)) nevers("not implemented");

    u a(ia), b(ib), m(imod);
    if (!((a >> 63) | (b >> 63)))
    {
        auto ab = a + b;
        return { ab % m, ab >= m };
    }

    nevers("not yet implemented using lifter");
}

static Integer mod_mul_rec(Integer a, Integer b, Integer m)
{
    if (b == 0) return 0;
    if (b == 1) return a % m;

    a %= m;
    b %= m;

    auto y = (a << 1);
    auto q = (b >> 1);
    y = mod_mul_rec(y, q, m);
    if (b == (q << 1)) return y;
    return (y + a) % m;
}

static Integer mod_mul_max(Integer a, Integer b, Integer m)
{
    if (a > b) return mod_mul_rec(a, b, m);
    return mod_mul_rec(b, a, m);
}

Integer modpow(Integer x, Integer p, Integer m)
{
    if ( p == 0 ) return 1;
    if ( p == 1 ) return modmul(x, 1, m);
    if ( p == 2 ) return modmul(x, x, m);

    auto y = modpow(x, p / 2, m);
    auto z = mod_mul_rec(y, y, m);

    if ( !(p % 2) ) return z;

    return mod_mul_max(x, z, m);
}

