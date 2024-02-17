#pragma once

#include <vector>

#include "err.h"
#include "integer.h"
#include "rns.h"

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
Poly rangeUpP(Poly a, Integer q); // moves negative values to q-range
Poly rangeDownP(Poly a, Integer q); // moves down values to q-range (if above q)
Poly rangeCenterP(Poly a, Integer q); // moves to (-q/2,q/2)-range

Poly neg(Poly a, Integer q);
Poly mulmod_simple(const Poly & a, const Poly & b, Integer q);
Poly mul_simple(const Poly & a, const Poly & b);
Poly mulmod_elwise(const Poly & a, const Poly & b, Integer q);
Poly mulmod_btrfly(const Poly & a, const Poly & b, Integer q);
//inline Poly mul(Poly a, Poly b, Integer q) { return mulmod_simple(a, b, q); }
inline Poly mul(Poly a, Poly b, Integer q)
{
    auto r = mulmod_btrfly(a, b, q);
    if (1)  // debug
    {
        auto s = mulmod_simple(a, b, q);
        if (r.v != s.v) never;
    }
    return r;
}
Poly add(Poly a, Poly b, Integer q);
Poly rescaleRound(const Poly & a, Integer idelta);
//Poly rescaleFloor(const Poly & a, Integer w);
Poly mul(Poly a, Integer b, Integer q);
Poly div(Poly a, Integer b, Integer q);
Poly scaleUp(const Poly & a, Integer Q, Integer P, Integer PQ);

} //poly

namespace mod
{

//inline Integer add(Integer a, Integer b, Integer mod)
//{
//
//inline Integer sub(Integer a, Integer b, Integer mod)
//{

} // mod

std::ostream & operator<<(std::ostream & os, const poly::Poly & p);


namespace poly
{

struct PolyRns
{
    ///static rns_ns::RnsMrs rns;
    using Tower = Poly;

    std::vector<Tower> towers;
    const rns_ns::Rns * rns_ptr;

    void operator+=(Integer x);
    void operator+=(rns_ns::RnsForm x);
    ///    int size() const { return (int)v.size(); }
    PolyRns(const rns_ns::Rns & r) : rns_ptr(&r) {}
    PolyRns(const rns_ns::Rns * p) : rns_ptr(p) {}
    //    explicit Poly(int n) : v(n) {}
    PolyRns(const rns_ns::Rns * r, int n, Integer i);
    Tower getCollapsed(bool forceBlend = false) const;

    int polysize() const { if (towers.empty())never; return towers[0].size(); }
    int ntows() const { return rns_ptr->size(); }
    rns_ns::RnsForm rnsForm(int i) const;
    void negInplace();
    bool match(const PolyRns & b) const;

    PolyRns rebase(const rns_ns::Rns & nr) const;
    PolyRns modswap(const rns_ns::Rns & nr) const;
    PolyRns shrink(const rns_ns::RnsShrinkRound & rshr) const;
};

PolyRns mul(const PolyRns & a, const PolyRns & b);
PolyRns mul(const PolyRns & a, rns_ns::RnsForm b);
//PolyRns mul(PolyRns a, Integer b);
//PolyRns mul(PolyRns a, const vint & vb);
PolyRns neg(PolyRns a);
PolyRns add(PolyRns a, PolyRns b);
PolyRns mul_simple(const PolyRns & a, const PolyRns & b);
PolyRns rescaleRoundRns(const PolyRns & a, Integer idelta);

} // poly

inline std::ostream & operator<<(std::ostream & os, const poly::PolyRns & p)
{
    return os << p.getCollapsed(true);
}
