#include <cmath>
#include <sstream>

#include <iostream> // debug
using std::cout;

#include "ckkshyb.h"
#include "err.h"
#include "mathut.h"

using poly::Poly;

ckks::EkHybP::EkHybP(int level, SkP sk, Param p, RndStream& rs)
{
    if (p.w == 0) nevers("Digit size is not set; assign size to 'w'");
    Poly s = sk.s;
    int n = (int)s.size();
    Integer in(n), iU(1);

    never;
    //ql = p.q(level);
    //P = p.w;

    //// ntt requires inversion n in PQ
    //for (int i = 0; i < 10000000; i++, ++P)
    //    if (math::gcdT(P, in) == iU)
    //        break;

    //auto PQl = P * ql;
    //const auto& q = PQl;

    //int dnum = poly::calc_dnum(p.w, q); // ql - doesnt work, error in paper

    ////Poly a = genPolyRq(sk.n, rs, q);
    ////Poly e = genPolyEr(sk.n, rs);

    //da.resize(dnum);
    //db.resize(dnum);
    //poly::Dpoly e(dnum);
    //for (int i = 0; i < dnum; i++)
    //{
    //    da[i] = genPolyRq(sk.n, rs, q);
    //    e[i] = genPolyEr(sk.n, rs);
    //    e[i] = rangeUp(e[i], q);
    //}

    //s = rangeUp(s, q);
    //auto s2 = mul(s, s, q);
    //auto ds2 = poly::PW(s2, p.w, q);

    //// b = -a*SK + e + P*SK*SK
    //// a = a
    //for (int i = 0; i < dnum; i++)
    //{
    //    auto x1 = neg(da[i], q);
    //    auto x2 = mul(x1, s, q);
    //    auto x3 = add(x2, e[i], q);
    //    Poly x5 = mul(ds2[i], P, q);
    //    db[i] = add(x3, x5, q);
    //}
}

string ckks::EkHybP::print() const
{
	never;
	return string();
}

poly::Dpoly poly::WDp(const poly::Poly& a, Integer w, Integer q)
{
	never;
	return Dpoly();
}

poly::Dpoly poly::PWp(const poly::Poly& a, Integer w, Integer q)
{
	never;
	return Dpoly();
}

int poly::calc_dnumP(Integer w, Integer q)
{
	never;
	return 0;
}

poly::Poly poly::dotP(const Dpoly& a, const Dpoly& b, Integer q)
{
	never;
	return poly::Poly();
}

ckks::CtxtP ckks::relinHybP(const Ctxt3P& c, const Param& p, const EkHybP& ek)
{
	never;
	//return CtxtP();
}

ckks::CtxtP ckks::mulHybP(const CtxtP& a, const CtxtP& b, const Param& p, const EkHybP& ek)
{
	never;
	//return CtxtP();
}
