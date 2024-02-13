#include "integer.h"

const int TESTNPRIMES = 3;

bool isPrime(Integer x);

Integer factor_ferma(Integer x);
Integer factor_root(Integer x);
Integer factor_pollard_pm1(Integer x);
Integer factor_pollard_rho(Integer x);

std::vector<Integer> factors(Integer x, int alg);
// alg 0:root, 1:ferma, 2:pollardPm1, 3:pollardRho
std::vector<Integer> factors(Integer x); // automatic selection

template <class T>
inline void operator+=(std::vector<T> & a, const std::vector<T> & b)
{ a.insert(a.end(), b.begin(), b.end()); }


// find max n such that (n+1)^2 > x
// newton root square integer
template <class T> T sqrtnT(T n)
{
    const T t2 {2}, t4 {4};

    if ( n < t2 ) return n;
    if ( n < t4 ) return (n + 1) / 2;

    T x = n, y = n;
    while ( y > t2 ) { x >>= 1; y >>= 2; }

    // x>1
    while ( x * x < n ) x += (x >> 1);

    for (;;)
    {
        auto z = (x + n / x) / 2;
        if ( z >= x ) break;
        x = z;
    }

    return x;
}

inline Integer sqrtn_newton(Integer x) { return sqrtnT<Integer>(x); }
inline Integer sqrtn(Integer x) { return sqrtn_newton(x); }
inline bool isSqrtn(Integer x) { auto s = sqrtn(x); return s * s == x; }
