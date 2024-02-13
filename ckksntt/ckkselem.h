#pragma once

#include <complex>
#include <vector>

#include "bigun.h"
#include "err.h"
#include "poly.h"

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
        int levels;
        std::vector<Integer> vql;
        ParamQx(Integer q, int lev, ParamEncode pe) : penc(pe), levels(lev)
        {
            vql.push_back(q);
            for (int i = 0; i < levels; i++) vql.push_back(q *= penc.idelta);
            forceNttValues();
        }
        Integer q0() const { return vql[0]; }
        Integer qL() const { return vql[levels]; }
        Integer q(int l) const { return vql[l]; }

    private:
        void forceNttValues();
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

poly::Poly encode(const Param & p, const std::vector<cx> & v);
std::vector<cx> decode(const Param & p, const poly::Poly & m);

class RndStream
{
        int r2 = 0;
        Integer rq { 0 };
        int er = 0;

    public:
        int getR2();
        Integer getRq(Integer q);
        Integer getEr();
};

struct Sk
{
    poly::Poly s;
    int n;
    Sk(RndStream & rs, int sz);
};

struct Ctxt
{
    using Poly = poly::Poly;
    int level;
    Poly c0, c1;
    Ctxt(int lev, const Poly & a0, const Poly & a1) : level(lev), c0(a0), c1(a1) {}
};

struct Ctxt3 : Ctxt
{
    poly::Poly c2;
    Ctxt3(const Ctxt & c01, const Poly & c2) : Ctxt(c01), c2(c2) {}
};

// FIXME RnsStream must be ref&
poly::Poly genPolyRq(int n, RndStream rs, Integer q);
poly::Poly genPolyEr(int n, RndStream rs);
poly::Poly genPolyR2(int n, RndStream rs);

struct Pk
{
    int n;
    poly::Poly p0, p1;
    Pk(Sk sk, Param p, RndStream rs);
};

Ctxt encrypt(Pk pk, poly::Poly m, Param p, RndStream rs);
Ctxt encrypt(Sk sk, poly::Poly m, Param p, RndStream rs);
poly::Poly decrypt(Sk sk, Ctxt c, Param p);
poly::Poly decrypt(Sk sk, Ctxt3 c, Param p);

Ctxt add(const Ctxt & a, const Ctxt & b, Param p);
Ctxt3 mul3(const Ctxt & a, const Ctxt & b, Param p);
Ctxt3 rescale(const Ctxt3 & c, Integer idelta);
Ctxt rescale(const Ctxt & c, Integer idelta);

struct EkExt
{
    Integer P;
    poly::Poly b, a;
    EkExt(Sk sk, Param p, RndStream rs);
};

Ctxt mulExt(const Ctxt & a, const Ctxt & b, const Param & p, const EkExt & ek);
Ctxt relinExt(const Ctxt3 & c, const Param & p, const EkExt & ek);

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
