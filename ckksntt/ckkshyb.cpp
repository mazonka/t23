#include <cmath>
#include <sstream>

#include <iostream> // debug
using std::cout;

#include "ckkshyb.h"
#include "err.h"
#include "mathut.h"

using poly::Poly;

ckks::EkHyb::EkHyb(int lev, Sk sk, Param p, RndStream & rs) : level(lev)
{
    if (p.w == 0) nevers("Digit size is not set; assign size to 'w'");
    Poly s = sk.s;
    int n = (int)s.size();
    Integer in(n), iU(1);

    ql = p.q(level);
    P = p.w;

    // ntt requires inversion n in PQ
    for (int i = 0; i < 10000000; i++, ++P)
        if (math::gcdT(P, in) == iU)
            break;

    auto PQl = P * ql;
    const auto & q = PQl;

    int dnum = poly::calc_dnum(p.w, q); // ql - doesnt work, error in paper

    //Poly a = genPolyRq(sk.n, rs, q);
    //Poly e = genPolyEr(sk.n, rs);

    da.resize(dnum);
    db.resize(dnum);
    poly::Dpoly e(dnum);
    for (int i = 0; i < dnum; i++)
    {
        da[i] = genPolyRq(sk.n, rs, q);
        e[i] = genPolyEr(sk.n, rs);
        e[i] = rangeUp(e[i], q);
    }

    s = rangeUp(s, q);
    auto s2 = mul(s, s, q);
    auto ds2 = poly::PW(s2, p.w, q);

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

string ckks::EkHyb::print() const
{
    std::ostringstream os;
    //int level;
    //Integer P, ql;
    //poly::Dpoly db, da;
    cout << "level = " << level << " \t| current level of relinearization\n";
    cout << "P = " << todbl(P) << " \t| extension; = approx w; must be coprime to n\n";
    cout << "ql = " << todbl(ql) << " \t| current Ql value\n";
    cout << "db =\n";
    for (auto x : db) cout << ' ' << x << '\n';
    cout << "da =\n";
    for (auto x : da) cout << ' ' << x << '\n';
    return os.str();
}

poly::Dpoly poly::WD(const poly::Poly & a, Integer w, Integer q)
{
    int dnum = calc_dnum(w, q);
    Dpoly r;
    Poly b = a;
    for (int i = 0; i < dnum; i++)
    {
        Poly c = rangeDown(b, w);
        r.push_back(c);
        b = rescaleFloor(b, w);
    }
    return r;
}

poly::Dpoly poly::PW(const poly::Poly & a, Integer w, Integer q)
{
    int dnum = calc_dnum(w, q);
    Dpoly r;
    Poly b = a;
    for (int i = 0; i < dnum; i++)
    {
        r.push_back(b);
        b = poly::mul(b, w, q);
    }
    return r;
}

int poly::calc_dnum(Integer w, Integer q)
{
    int r = 0;
    while (q != Integer(0)) { q /= w; ++r; }
    return r;
}

poly::Poly poly::dot(const Dpoly & a, const Dpoly & b, Integer q)
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

ckks::Ctxt ckks::relinHyb(const Ctxt3 & c, const Param & par, const EkHyb & ek)
{
    // reference: pages 6,7,8 Polyakov 2021 Approximate...
    Ctxt r = c; // slice object

    Integer q = par.q(c.level);
    if (q != ek.ql) nevers("wrong relin Q level");

    auto w = par.w;
    auto P = ek.P;
    Integer pq = P * q;

    //auto d2 = scaleUp(c.c2, q, P, pq);
    auto d2 = rangeUp(c.c2, q);

    poly::Dpoly wd2 = poly::WD(d2, w, pq);
    auto d2eka = poly::dot(wd2, ek.da, pq);
    auto d2ekb = poly::dot(wd2, ek.db, pq);
    auto pa = div(d2eka, ek.P, pq);
    auto pb = div(d2ekb, ek.P, pq);
    r.c0 = add(r.c0, pb, q);
    r.c1 = add(r.c1, pa, q);

    return r;
}

ckks::Ctxt ckks::mulHyb(const Ctxt & a, const Ctxt & b, const Param & p, const EkHyb & ek)
{
    Ctxt3 c3 = mul3(a, b, p);
    Ctxt c2 = relinHyb(c3, p, ek);
    Ctxt c2sc = rescale(c2, p.penc.idelta);
    return c2sc;
}
