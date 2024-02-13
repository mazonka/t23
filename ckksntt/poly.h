#pragma once

#include <vector>

#include "bigun.h"
#include "err.h"

const bool DB = true;

namespace poly
{

struct Poly
{
    std::vector<Integer> v;
    void operator+=(Integer x) { v.push_back(x); }
    int size() const { return (int)v.size(); }
    Poly() {}
    explicit Poly(int n) : v(n) {}
    Poly(int n, Integer i) : v(n, i) {}
};

// Poly operations
Poly rangeUp(Poly a, Integer q); // moves negative values to q-range
Poly rangeDown(Poly a, Integer q); // moves down values to q-range (if above q)
Poly rangeCenter(Poly a, Integer q); // moves to (-q/2,q/2)-range

Poly neg(Poly a, Integer q);
Poly mulmod_simple(const Poly & a, const Poly & b, Integer q);
Poly mul_simple(const Poly & a, const Poly & b);
Poly mulmod_elwise(const Poly & a, const Poly & b, Integer q);
Poly mulmod_btrfly(const Poly & a, const Poly & b, Integer q);
//inline Poly mul(Poly a, Poly b, Integer q) { return mulmod_simple(a, b, q); }
inline Poly mul(Poly a, Poly b, Integer q)
{
    auto r = mulmod_btrfly(a, b, q);
    if (DB)
    {
        auto s = mulmod_simple(a, b, q);
        if (r.v != s.v) never;
    }
    return r;
}
Poly add(Poly a, Poly b, Integer q);
Poly rescaleRound(const Poly & a, Integer idelta);
Poly rescaleFloor(const Poly & a, Integer w);
Poly mul(Poly a, Integer b, Integer q);
Poly div(Poly a, Integer b, Integer q);
Poly scaleUp(const Poly & a, Integer Q, Integer P, Integer PQ);

} //poly

namespace mod
{

inline Integer add(Integer a, Integer b, Integer mod)
{
    auto r = a + b;
    if (r >= mod) r -= mod;
    if (DB) if (r >= mod) nevers("mod overflow");
    return r;
}

inline Integer sub(Integer a, Integer b, Integer mod)
{
    if (a >= b) return a - b;
    auto r = a + mod - b;
    if (DB) if (r < 0) nevers("mod overflow");
    return r;
}

} // mod

std::ostream & operator<<(std::ostream & os, const poly::Poly & p);
template <class T>
inline std::ostream & operator<<(std::ostream & os, const std::vector<T> & v)
{
    for (auto x : v) os << ' ' << x;
    return os;
}

