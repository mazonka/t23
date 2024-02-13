#include <algorithm>

#include "ccrun.h"
//#include "ccrut.h"

using Integer = int;
using vint = std::vector<Integer>;

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
    if (gcd * gcd != 1 ) nevers("no inverse (ull)");
    while (s < 0) s += mod;
    s /= gcd;
    return s % mod;
}



void cmain_v3()
{
    int a, m = 19, b = 256;
    int m1b = invert(m, b);
    int m2 = m / 2;

    for ( int i = 0; i < m; i++ )
    {
        a = i;
        int e1 = int(a * b * 1.0 / m + 0.5);
        int e2 = (m1b * (m2 + b - ( (a * b + m2) % m )) % b) % b;
        if ( e1 != e2 ) never;
        cout << i << '\t' << e1 << '\t' << e2 << '\n';
    }

}

void cmain_v4()
{
    int a, m = 19, b = 256;
    int m1mb = b - invert(m, b);

    for ( int i = 0; i < m; i++ )
    {
        a = i;
        int e1 = int(a * b * 1.0 / m);
        int e2 = ( ( (a * b) % m ) * m1mb ) % b;
        if ( e1 != e2 ) never;
        cout << i << '\t' << e1 << '\t' << e2 << '\n';
    }

}

void cmain()
{
    cmain_v3();
    cmain_v4();
}
