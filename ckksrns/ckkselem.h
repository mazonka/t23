#pragma once

#include <complex>
#include <vector>

#include "integer.h"
#include "err.h"
#include "poly.h"
#include "rns.h"

using cx = std::complex<double>;

namespace ckks
{

struct ParamEncode
{
    int n;
    Integer idelta;
    double ddelta;

    std::vector<cx> vxi;

    void init();
    ParamEncode(int an, Integer delta) : n(an), idelta(delta), ddelta(todbl(delta)) { init(); }

    cx mxA(int i, int j) const { return vxi[(j * (2 * i + 1)) % (2 * n)]; }
    cx mxAm1(int i, int j) const { return conj(mxA(j, i)); }
};

struct ParamQx // parameters for ModUp ModDown
{
        ParamEncode penc;
        int levels; // level+1==vqs.size()
        std::vector<Integer> vqs;
        ParamQx(Integer q, int lev, ParamEncode pe) : penc(pe), levels(lev)
        {
            vqs = findNttCoprimes(q);
            std::reverse(vqs.begin(), vqs.end());
            // reverse is done to move the best (closest) to the
            // end, so error on rescaling is minimal
        }
        Integer q0() const { return vqs[0]; }
        ///Integer qL() const { return vql[levels]; }
        ///Integer q(int l) const { return vql[l]; }
        Integer qL_() const;
        Integer q_(int l) const;
        rns_ns::RnsForm qLrns(const rns_ns::Rns & rns) const;

    private:
        std::vector<Integer> findNttCoprimes(Integer q);
};

struct ParamHyb
{
    Integer w { 0 }; // digit; e.g. 64 or 256
};

struct Param : ParamQx, ParamHyb
{
    Param(int n, Integer q0, Integer delta, int levs)
        : ParamQx(q0, levs, ParamEncode(n, delta)) {}
    string print() const;
};

poly::Poly encodeP(const Param & p, const std::vector<cx> & v);
poly::PolyRns encodeR(const Param & p, const std::vector<cx> & v, const rns_ns::Rns & rns);
std::vector<cx> decodeP(const Param & p, const poly::Poly & m);
std::vector<cx> decodeR(const Param & p, const poly::PolyRns & m, const rns_ns::Rns & rns);

class RndStream
{
        int r2 = 0;
        Integer rq { 0 };
        int er = 0;

    public:
        int getR2();
        Integer getRqSeed();
        Integer getEr();

        RndStream() {};
        RndStream(const RndStream &) = delete;
        RndStream & operator=(const RndStream &) = default;
};

struct SkP
{
    poly::Poly s;
    int n;
    SkP(RndStream & rs, int sz);
};

struct SkR
{
    int n;
    poly::PolyRns s;
    SkR(RndStream & rs, int sz, const rns_ns::Rns & r);
};

struct CtxtP
{
    using Poly = poly::Poly;
    int level;
    Poly c0, c1;
    CtxtP(int lev, const Poly & a0, const Poly & a1) : level(lev), c0(a0), c1(a1) {}
};

struct CtxtR
{
    using PolyR = poly::PolyRns;
    int level;
    PolyR c0, c1;
    CtxtR(int lev, const PolyR & a0, const PolyR & a1) : level(lev), c0(a0), c1(a1) {}
};

struct Ctxt3P : protected CtxtP
{
    poly::Poly c2;
    Ctxt3P(const CtxtP & c01, const Poly & c2) : CtxtP(c01), c2(c2) {}

    CtxtP slice() const { return *this;  }
    using CtxtP::level;
    using CtxtP::c0, CtxtP::c1;
};

struct Ctxt3R : protected CtxtR
{
    using PolyR = poly::PolyRns;
    PolyR c2;
    Ctxt3R(const CtxtR & c01, const PolyR & c2) : CtxtR(c01), c2(c2) {}

    CtxtR slice() const { return *this; }
    using CtxtR::level;
    using CtxtR::c0, CtxtR::c1;
};

// FIXME RnsStream must be ref&
poly::Poly genPolyRqP(int n, RndStream & rs, Integer q);
poly::Poly genPolyErP(int n, RndStream & rs);
poly::Poly genPolyR2P(int n, RndStream & rs);

poly::PolyRns genPolyRqR(int n, RndStream & rs, rns_ns::RnsForm q, const rns_ns::Rns & r);
poly::PolyRns genPolyErR(int n, RndStream & rs, const rns_ns::Rns & r);
poly::PolyRns genPolyR2R(int n, RndStream & rs, const rns_ns::Rns & r);

rns_ns::RnsForm getRqRns(RndStream & rs, const rns_ns::RnsForm & q);
Integer getRq(RndStream & rs, Integer q);


struct PkP
{
    int n;
    poly::Poly p0, p1;
    PkP(SkP sk, Param p, RndStream & rs);
};

struct PkR
{
    int n;
    poly::PolyRns p0, p1;
    PkR(SkR sk, Param p, RndStream & rs);
};

CtxtP encryptP(SkP sk, poly::Poly m, Param p, RndStream & rs);
CtxtP encryptP(PkP pk, poly::Poly m, Param p, RndStream & rs);
CtxtR encryptR(SkR sk, poly::PolyRns m, Param p, RndStream & rs);
CtxtR encryptR(PkR pk, poly::PolyRns m, Param p, RndStream & rs);
poly::Poly decryptP(SkP sk, CtxtP c, Param p);
poly::Poly decryptP3(SkP sk, Ctxt3P c, Param p);
poly::PolyRns decryptR(SkR sk, CtxtR c, Param p);
poly::PolyRns decryptR3(SkR sk, Ctxt3R c, Param p);
//poly::Poly decrypt(Sk sk, Ctxt3 c, Param p);


CtxtP add(const CtxtP & a, const CtxtP & b, Param p);
CtxtR add(const CtxtR & a, const CtxtR & b);
Ctxt3P mul3(const CtxtP & a, const CtxtP & b, Param p);
Ctxt3R mul3(const CtxtR & a, const CtxtR & b);
Ctxt3P rescale(const Ctxt3P & c, Integer idelta, Param par);
Ctxt3R rescale(const Ctxt3R & c, Integer idelta);
CtxtP rescale(const CtxtP & c, Integer idelta, Param par);
CtxtP rescaleLevel(const CtxtP & c, Param par);
CtxtR rescale(const CtxtR & c, Integer idelta, Param par);
CtxtR rescaleLevel(const CtxtR & c, const rns_ns::RnsShrinkRound & dat);

struct EkExtP
{
    Integer P;
    poly::Poly b, a;
    EkExtP(SkP sk, Param p, RndStream & rs, Integer newp = 0);
};

struct EkExtR
{
    Integer P;
    poly::PolyRns b, a; // are in rext
    rns_ns::RnsShrinkRound rshrink;
    EkExtR(SkR sk, Param p, RndStream & rs,
           rns_ns::Rns & rext, rns_ns::RnsShrinkRound rshrink);
};

CtxtP relinExt(const Ctxt3P & c, const Param & p, const EkExtP & ek);
CtxtR relinExt(const Ctxt3R & c, const Param & p, const EkExtR & ek);
CtxtP mulExtP(const CtxtP& a, const CtxtP& b, const Param& p, const EkExtP& ek);
CtxtR mulExtR(const CtxtR& a, const CtxtR& b, const Param& p, const EkExtR& ek, const rns_ns::RnsShrinkRound& datQ);

} // ckks


double roundd(double prec, double x);
inline cx roundc(double p, cx a) { return cx { roundd(p, a.real()), roundd(p, a.imag()) }; }

template <class T> inline T roundx(double p, T x) { never; }
template <> inline double roundx<double>(double p, double x) { return roundd(p, x); }
template <> inline cx roundx<cx>(double p, cx x) { return roundc(p, x); }

template <class T>
inline std::vector<T> roundv(double p, const std::vector<T> & v)
{
    auto r = v;
    for (auto & x : r) x = roundx(p, x);
    return r;
}


std::vector<cx> operator*(const std::vector<cx> & v, double x);
