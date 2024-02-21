#include <cmath>
#include <sstream>

#include <iostream> // debug
using std::cout;

#include "ckkshyb.h"
#include "err.h"
#include "mathut.h"
#include "egcd.inc"

using poly::Poly;

ckks::EkHybP::EkHybP(int level, SkP sk, Param p, RndStream & rs)
{
    if (p.w == 0) nevers("Digit size is not set; assign size to 'w'");
    Poly s = sk.s;
    int n = (int)s.size();
    Integer in(n), iU(1);

    ql = p.q_(level);
    P = p.w;

    // ntt requires inversion n in PQ
    for (int i = 0; i < 10000000; i++, ++P)
        if (gcdT(P, in) == iU)
            break;

    auto PQl = P * ql;
    const auto& q = PQl;

    int dnum = poly::calc_dnumP(p.w, q); // ql - doesnt work, error in paper

    //Poly a = genPolyRq(sk.n, rs, q);
    //Poly e = genPolyEr(sk.n, rs);

    da.resize(dnum);
    db.resize(dnum);
    poly::Dpoly e(dnum);
    for (int i = 0; i < dnum; i++)
    {
        da[i] = genPolyRqP(sk.n, rs, q);
        e[i] = genPolyErP(sk.n, rs);
        e[i] = rangeUpP(e[i], q);
    }

    s = rangeUpP(s, q);
    auto s2 = mul(s, s, q);
    auto ds2 = poly::PWp(s2, p.w, q);

    // b = -a*SK + e + P*SK*SK
    // a = a
    for (int i = 0; i < dnum; i++)
    {
        auto x1 = neg(da[i], q);
        auto x2 = mul(x1, s, q);
        auto x3 = add(x2, e[i], q);
        Poly x5 = mul(ds2[i], P, q);
        db[i] = add(x3, x5, q);
    }
}

string ckks::EkHybP::print() const
{
    never;
    return string();
}

poly::Dpoly poly::WDp(const poly::Poly & a, Integer w, Integer q)
{
    int dnum = calc_dnumP(w, q);
    Dpoly r;
    Poly b = a;
    for (int i = 0; i < dnum; i++)
    {
        Poly c = rangeDownP(b, w);
        r.push_back(c);
        b = rescaleFloor(b, w);
    }
    return r;
}

poly::Dpoly poly::PWp(const poly::Poly & a, Integer w, Integer q)
{
    int dnum = calc_dnumP(w, q);
    Dpoly r;
    Poly b = a;
    for (int i = 0; i < dnum; i++)
    {
        r.push_back(b);
        b = poly::mul(b, w, q);
    }
    return r;
}

int poly::calc_dnumP(Integer w, Integer q)
{
    int r = 0;
    while (q != Integer(0)) { q /= w; ++r; }
    return r;
}

poly::Poly poly::dotP(const Dpoly & a, const Dpoly & b, Integer q)
{
    int dnum = (int)a.size();
    auto n = a[0].size();
    Poly r(n, Integer(0));
    for (int i = 0; i < dnum; ++i)
    {
        Poly s = poly::mul(a[i], b[i], q);
        r = poly::add(r, s, q);
    }
    return r;
}

ckks::CtxtP ckks::relinHybP(const Ctxt3P & c, const Param & par, const EkHybP & ek)
{
    // reference: pages 6,7,8 Polyakov 2021 Approximate...
    CtxtP r = c.slice(); // slice object

    Integer q = par.q_(c.level);
    if (q != ek.ql) nevers("wrong relin Q level");

    auto w = par.w;
    auto P = ek.P;
    Integer pq = P * q;

    //auto d2 = scaleUp(c.c2, q, P, pq);
    auto d2 = rangeUpP(c.c2, q);

    poly::Dpoly wd2 = poly::WDp(d2, w, pq);
    auto d2eka = poly::dotP(wd2, ek.da, pq);
    auto d2ekb = poly::dotP(wd2, ek.db, pq);
    auto pa = div(d2eka, ek.P, pq);
    auto pb = div(d2ekb, ek.P, pq);
    r.c0 = add(r.c0, pb, q);
    r.c1 = add(r.c1, pa, q);

    return r;
}

ckks::CtxtP ckks::mulHybP(const CtxtP & a, const CtxtP & b, const Param & p, const EkHybP & ek)
{
    Ctxt3P c3 = mul3(a, b, p);
    CtxtP c2 = relinHybP(c3, p, ek);
    CtxtP c2sc = rescaleLevel(c2, p);
    return c2sc;
}
