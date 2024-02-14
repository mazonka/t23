///#include <iostream> /// debug
///using std::cout;


#include "poly.h"
#include "ntt.h"


poly::Poly poly::rangeUp(Poly a, Integer q)
{
    Poly r = a;
    for (auto & x : r.v) while ( x < 0 ) x += q;
    return r;
}

poly::Poly poly::rangeDown(Poly a, Integer q)
{
    Poly r = a;
    for (auto & x : r.v) x = x % q;
    return r;
}

poly::Poly poly::rangeCenter(Poly a, Integer q)
{
    Poly r = a;
    auto q2 = q / 2;
    for (auto & x : r.v) if ( x > q2 ) x -= q;
    return r;
}

poly::Poly poly::neg(Poly a, Integer q)
{
    Poly r;
    for (auto x : a.v) r += q - x;
    return r;
}

poly::Poly poly::mulmod_simple(const Poly & a, const Poly & b, Integer q)
{
    int cntr = 0;

    auto n = b.size();  if (a.size() != n) never;
    Poly r(n * 2, Integer(0));

    auto assertRange = [&](Integer x)
    {
        //cout << (++cntr) << " " << x << '\n';
        if (x < Integer(0)) never;
        if (x >= q) never;
    };

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            auto & rij = r.v[i + j];
            ///cout << (++cntr) << " Rij a b : " << rij <<' '<< a.v[i]<<' '<< b.v[j] << '\n';
            rij += mod::mul(a.v[i], b.v[j], q);
            ///cout << (++cntr) << " Rij : " << rij << '\n';
            if (rij >= q) rij -= q;
            //assertRange(rij);
        }

    for (int i = 0; i < n; i++)
    {
        auto & d = r.v[i];
        auto & u = r.v[i + n];
        if (d < u) d += q;
        d -= u;
        //assertRange(d);
    }
    r.v.resize(n);

    ///cout << "c " << r << '\n';
    return r;
}

poly::Poly poly::mul_simple(const Poly & a, const Poly & b)
{
    auto n = b.size();  if (a.size() != n) never;
    Poly r(n * 2, Integer(0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            auto & rij = r.v[i + j];
            rij += a.v[i] * b.v[j];
        }

    for (int i = 0; i < n; i++)
    {
        auto & d = r.v[i];
        auto & u = r.v[i + n];
        d -= u;
    }
    r.v.resize(n);

    return r;
}

poly::Poly poly::mulmod_elwise(const Poly & a, const Poly & b, Integer q)
{
    auto sz = b.size();  if (a.size() != sz) never;
    Poly r(sz);
    for (int i = 0; i < sz; i++) r.v[i] = mod::mul(a.v[i], b.v[i], q);
    return r;
}

poly::Poly poly::mulmod_btrfly(const Poly & a, const Poly & b, Integer q)
{
    if (ntt::disabled) return mulmod_simple(a, b, q);

    using namespace ntt;
    auto an = nttBfly(a, q);
    auto bn = nttBfly(b, q);
    auto cn = mulmod_elwise(an, bn, q);
    auto c = ittBfly(cn, q);
    return c;
}


poly::Poly poly::add(Poly a, Poly b, Integer q)
{
    auto sz = b.size();  if (a.size() != sz) never;

    Poly r(sz);
    for (int i = 0; i < sz; i++)
    {
        r.v[i] = a.v[i] + b.v[i];
        if (r.v[i] >= q) r.v[i] -= q;
        if (1) if (r.v[i] >= q) nevers("input overflow");
    }
    return r;
}

poly::Poly poly::rescaleRound(const Poly & a, Integer idelta)
{
    auto d2 = idelta / 2;
    Poly r(a);
    for (auto & x : r.v) x = (x + d2) / idelta;
    return r;
}

poly::Poly poly::rescaleFloor(const Poly & a, Integer w)
{
    Poly r(a);
    for (auto & x : r.v) x /= w;
    return r;
}

poly::Poly poly::mul(Poly a, Integer b, Integer q)
{
    auto sz = a.size();

    Poly r(sz);
    for (int i = 0; i < sz; i++)
        r.v[i] = mod::mul(a.v[i], b, q);

    return r;
}

poly::Poly poly::div(Poly a, Integer b, Integer q)
{
    a = rescaleRound(a, b);
    a = rangeDown(a, q);
    return a;
}

poly::Poly poly::scaleUp(const Poly & a, Integer Q, Integer P, Integer PQ)
{
    Poly r = rangeCenter(a, Q);
    for (auto & x : r.v) x = x * P;
    r = rangeUp(r, PQ);
    return r;
}

std::ostream & operator<<(std::ostream & os, const poly::Poly & p)
{
    os << "{";
    for (auto x : p.v) os << ' ' << x;
    os << " }";
    return os;
}

