#include <cmath>
#include <sstream>
#include <algorithm>

#include <iostream> // debug
using std::cout;

#include "ckkshyb.h"
#include "err.h"
#include "mathut.h"
#include "egcd.inc"

using poly::Poly;
using poly::PolyRns;

ckks::EkHybP::EkHybP(int lev, SkP sk, Param p, RndStream & rs) : level(lev)
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
    const auto & q = PQl;

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

Integer ckks::EkHybR::findExtDigit(const vint & qs, int n)
{
    Integer in(n), iU(1);

    Integer P = *std::max_element(qs.begin(), qs.end());
    ++P;

    // find P for extension
    for (int i = 0; i < 10000000; i++, ++P)
    {
        if (gcdT(P, in) != iU) continue; // ntt requires inversion n in PQ
        if (!isPrime(P)) continue;
        break;
    }
    return P;
}

ckks::EkHybR::EkHybR(SkR sk, Param p, RndStream & rs,
                     rns_ns::Rns & rext, rns_ns::RnsShrinkRound rshr)
///: level(lev)
    : da(rext), db(rext), rshrink(rshr)
{
    ///if (p.w == 0) nevers("Digit size is not set; assign size to 'w'");
    PolyRns s = sk.s;
    ///int n = s.polysize();
    ///Integer in(n), iU(1);

    ///ql = p.q_(level);
    ///P = p.w;

    ///vint qs = p.vqs;
    ///qs.resize(level+1);

    ///P = *std::max_element(qs.begin(),qs.end());
    ///++P;

    // find P for extension
    //for (int i = 0; i < 10000000; i++, ++P)
    //{
    //    if (gcdT(P, in) != iU) continue; // ntt requires inversion n in PQ
    //    if (!isPrime(P)) continue;
    //    break;
    //}

    ///qs.push_back(P);

    //auto PQl = P * ql;
    //const auto& q = PQl;

    //int dnum = poly::calc_dnumP(p.w, q); // ql - doesnt work, error in paper
    int dnum = rext.size();
    auto qs = rext.getQs();

    //poly::PolyRns e(rext);
    //da.towers.resize(dnum);
    db.towers.resize(dnum);
    //e.towers.resize(dnum);

    rns_ns::RnsForm qrf(rext, 0); // all range
    da = genPolyRqR(sk.n, rs, qrf, rext);
    PolyRns e = genPolyErR(sk.n, rs, rext);

    //for (int i = 0; i < dnum; i++)
    //{
    //    da.towers[i] = genPolyRqP(sk.n, rs, qs[i]);
    //    e.towers[i] = genPolyErP(sk.n, rs);
    //    e.towers[i] = rangeUpP(e.towers[i], qs[i]);
    //}


    //auto s2 = mul(s, s);
    //auto ds2 = poly::PWr(s2);

    //// b = -a*SK + e + P*SK*SK
    //// a = a
    //for (int i = 0; i < dnum; i++)
    //{
    //    auto q = qs[i];
    //    auto x1 = neg(da.towers[i], q);
    //    auto x2 = mul(x1, s, q);
    //    auto x3 = add(x2, e[i], q);
    //    Poly x5 = mul(ds2[i], P, q);
    //    db[i] = add(x3, x5, q);
    //}

    auto se = s.modswap(rext);
    auto s2 = mul(se, se);
    auto ds2 = poly::PWr(s2);

    auto PinPQ = rshrink.PinQ.rebaseAdd(rext);

    // b = -a*SK + e + P*SK*SK
    // a = a
    auto x1 = neg(da);
    auto x2 = mul(x1, se);
    auto x3 = add(x2, e);
    ///auto x4 = mul(se, se);
    PolyRns x5 = mul(ds2, PinPQ);
    db = add(x3, x5);

    cout << "AAA " << __func__ << " PQ=" << rext.dynrange_() << " s=" << se << " a=" << da << " e=" << e << " b=" << db << '\n';
    cout << "AAA2 " << "x1x2x3(s2,ds2)x5 " << x1 << x2 << x3 << s2 << ds2 << x5 << '\n';

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

PolyRns poly::WDr(const PolyRns & a)
{
    const rns_ns::Rns * rns = a.rns_ptr;
    vint qs = rns->getQs();
    //vint ms = rns->getMs();
    vint us = rns->getUs();

    int ntows = a.ntows();
    int sz = a.polysize();

    PolyRns r(rns);

    for (int i = 0; i < ntows; i++)
    {
        const PolyRns::Tower & t = a.towers[i];
        PolyRns::Tower u;
        for (int j = 0; j < sz; j++)
        {
            Integer a = modmul(t.v[j], us[i], qs[i]);
            u.v.push_back(a);
        }
        r.towers.push_back(u);
    }

    return r;
}

PolyRns poly::PWr(const PolyRns & a)
{
    const rns_ns::Rns * rns = a.rns_ptr;
    vint qs = rns->getQs();
    vint ms = rns->getMs();
    //vint us = rns->getUs();

    int ntows = a.ntows();
    int sz = a.polysize();

    PolyRns r(rns);

    for (int i = 0; i < ntows; i++)
    {
        const PolyRns::Tower & t = a.towers[i];
        PolyRns::Tower u;
        for (int j = 0; j < sz; j++)
        {
            Integer a = modmul(t.v[j], ms[i], qs[i]);
            u.v.push_back(a);
        }
        r.towers.push_back(u);
    }

    return r;
}

PolyRns poly::dotR(const PolyRns & a, const PolyRns & b)
{
    const rns_ns::Rns * rns = a.rns_ptr;
    vint qs = rns->getQs();
    int ntows = a.ntows();
    if (ntows != b.ntows()) never;

    auto n = a.polysize();

    PolyRns r(rns);
    for (int i = 0; i < ntows; ++i)
    {
        Poly s = poly::mul(a.towers[i], b.towers[i], qs[i]);
        ///r = poly::add(r, s, q);
        r.towers.push_back(s);
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

ckks::CtxtR ckks::relinHybR(const Ctxt3R & c, const Param & p, const EkHybR & ek)
{
    CtxtR r = c.slice(); // slice object

    ///Integer q = par.q_(c.level);
    ///if (q != ek.ql) nevers("wrong relin Q level");

    ///auto w = par.w;
    ///auto P = ek.P;
    ///Integer pq = P * q;

    //auto d2 = scaleUp(c.c2, q, P, pq);
    ///auto d2 = rangeUpP(c.c2, q);
    const auto & rext = *ek.da.rns_ptr;
    auto d2 = c.c2.rebase(rext);

    //auto d2eka = mul(d2, ek.a);
    //auto d2ekb = mul(d2, ek.b);
    PolyRns wd2 = poly::WDr(d2);
    auto d2eka = poly::dotR(wd2, ek.da);
    auto d2ekb = poly::dotR(wd2, ek.db);

    ///auto pa = div(d2eka, ek.P, pq);
    ///auto pb = div(d2ekb, ek.P, pq);
    auto pa = d2eka.shrink(ek.rshrink);
    auto pb = d2ekb.shrink(ek.rshrink);
    r.c0 = add(r.c0, pb);
    r.c1 = add(r.c1, pa);
    cout << "AAA " << __func__ << " d2=" << d2
         << " ek:a:b=" << ek.da << ek.db
         << " d2ek=" << d2eka << d2ekb << '\n';
    return r;

}

ckks::CtxtP ckks::mulHybP(const CtxtP & a, const CtxtP & b, const Param & p, const EkHybP & ek)
{
    Ctxt3P c3 = mul3(a, b, p);
    CtxtP c2 = relinHybP(c3, p, ek);
    CtxtP c2sc = rescaleLevelP(c2, p);
    return c2sc;
}

ckks::CtxtR ckks::mulHybR(const CtxtR & a, const CtxtR & b, const Param & p, const EkHybR & ek, const rns_ns::RnsShrinkRound & datQ)
{
    Ctxt3R c3 = mul3(a, b);
    CtxtR c2 = relinHybR(c3, p, ek);
    CtxtR c2sc = rescaleLevelR(c2, datQ);
    return c2sc;
}


