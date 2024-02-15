#include <cmath>
#include <sstream>

#include <iostream> // debug
using std::cout;

#include "ckkselem.h"
#include "err.h"
#include "ntt.h"

using poly::Poly;

poly::Poly ckks::encode(const Param & p, const std::vector<cx> & v)
{
    std::vector<Integer> m;
    auto vl = v;
    int sz = (int)v.size();
    auto n = p.penc.n;
    if (2 * sz != n) never;

    for (int i = 0; i < sz; i++)
        vl.push_back(conj(v[sz - 1 - i]));

    std::vector<cx> m1(n, cx {});
    for (int i = 0; i < n; i++)
    {
        cx sum = 0;
        for (int j = 0; j < n; j++)
        {
            auto elA = p.penc.mxAm1(i, j);
            sum += elA * vl[j];
        }
        m1[i] = sum;
    }

    //cout << "AAA1 " << roundv(1e-3,m1) << '\n';

    auto m2 = m1 * (1.0 / n);

    //cout << "AAA2 " << roundv(1e-3, m2) << '\n';

    std::vector<double> m3;
    for (auto x : m2) m3.push_back(x.real());

    //cout << "AAA3 " << roundv(1e-3, m3) << '\n';

    std::vector<double> m4;
    for (auto x : m3) m4.push_back(x * p.penc.ddelta);

    //cout << "AAA4 " << roundv(1e-3, m4) << '\n';

    //std::vector<Integer> m5;
    Poly m5;
    for (auto x : m4) m5 += Integer((signed long long)roundd(1.0, x));

    //cout << "AAA5 " << m5 << '\n';

    return { m5 };
}

std::vector<cx> ckks::decode(const Param & p, const Poly & m)
{
    auto n = p.penc.n;

    std::vector<double> m1;
    for (auto x : m.v) m1.push_back(todbl(x));

    std::vector<double> m2 = m1;
    for (auto & x : m2) x /= p.penc.ddelta;

    int sz = n / 2;
    std::vector<cx> r(sz);

    for (int i = 0; i < sz; i++)
    {
        cx sum = 0;
        for (int j = 0; j < n; j++)
        {
            int idx = (j * (2 * i + 1)) % (2 * n);
            sum += p.penc.vxi[idx] * m2[j];
        }
        r[i] = sum;
    }

    return r;
}

ckks::Ctxt ckks::encrypt(Sk sk, Poly m, Param p, RndStream & rs)
{
    auto q = p.qL();

    Poly s = sk.s;
    Poly a = genPolyRq(sk.n, rs, q);
    Poly e = genPolyEr(sk.n, rs);

    s = rangeUp(s, q);
    e = rangeUp(e, q);
    auto mR = rangeUp(m, q);

    // c0 = -a*SK + e + m
    // c1 = a
    auto x1 = neg(a, q);
    auto x2 = mul(x1, s, q);
    auto x3 = add(x2, e, q);
    auto x4 = add(x3, mR, q);

    Ctxt r(p.levels, x4, a);
    return r;
}

ckks::Ctxt ckks::encrypt(Pk pk, Poly m, Param p, RndStream & rs)
{
    // c0 = PK0*u+e1+M
    // c1 = PK1*u+e2
    auto q = p.qL();
    Poly e0 = genPolyEr(pk.n, rs);
    Poly e1 = genPolyEr(pk.n, rs);
    Poly u = genPolyR2(pk.n, rs);

    u = rangeUp(u, q);
    e0 = rangeUp(e0, q);
    e1 = rangeUp(e1, q);
    auto mR = rangeUp(m, q);

    auto x1 = mul(pk.p0, u, q);
    auto x2 = add(x1, e0, q);
    auto x3 = add(x2, mR, q);
    auto x4 = mul(pk.p1, u, q);
    auto x5 = add(x4, e1, q);

    Ctxt r(p.levels, x3, x5);
    return r;
}

poly::Poly ckks::decrypt(Sk sk, Ctxt c, Param p)
{
    auto q = p.q(c.level);

    Poly s = sk.s;

    s = rangeUp(s, q);
    auto c0 = rangeDown(c.c0, q);
    auto c1 = rangeDown(c.c1, q);

    // m = c0 + c1*SK mod q0
    auto x1 = mul(c1, s, q);
    auto x2 = add(c0, x1, q);

    auto m = rangeCenter(x2, q);

    return m;
}

poly::Poly ckks::decrypt(Sk sk, Ctxt3 c, Param p)
{
    auto q = p.q(c.level);

    Poly s = sk.s;

    s = rangeUp(s, q);
    auto c0 = rangeDown(c.c0, q);
    auto c1 = rangeDown(c.c1, q);
    auto c2 = rangeDown(c.c2, q);

    // m = c0 + c1*SK + c2*SK^2 mod q0
    auto x1 = mul(c1, s, q);
    auto x2 = add(c0, x1, q);
    auto x3 = mul(c2, s, q);
    auto x4 = mul(x3, s, q);
    auto x5 = add(x2, x4, q);

    auto m = rangeCenter(x5, q);

    return m;
}

ckks::Ctxt ckks::add(const Ctxt & a, const Ctxt & b, Param p)
{
    int lv = a.level;
    if (b.level != lv) never;
    auto q = p.q(lv);
    Poly c0 = poly::add(a.c0, b.c0, q);
    Poly c1 = poly::add(a.c1, b.c1, q);
    return Ctxt(lv, c0, c1);
}

ckks::Ctxt3 ckks::mul3(const Ctxt & a, const Ctxt & b, Param p)
{
    int lv = a.level;
    if (b.level != lv) never;
    auto q = p.q(lv);

    Poly c0 = mul(a.c0, b.c0, q);
    Poly c1x1 = mul(a.c0, b.c1, q);
    Poly c1x2 = mul(a.c1, b.c0, q);
    Poly c2 = mul(a.c1, b.c1, q);
    Poly c1 = add(c1x1, c1x2, q);
    Ctxt c01(lv, c0, c1);
    return Ctxt3(c01, c2);
}

ckks::Ctxt3 ckks::rescale(const Ctxt3 & c, Integer idelta)
{
    int lv = c.level - 1;
    if (lv < 0) nevers("negative level");
    Poly c0 = rescaleRound(c.c0, idelta);
    Poly c1 = rescaleRound(c.c1, idelta);
    Poly c2 = rescaleRound(c.c2, idelta);
    return Ctxt3(Ctxt(lv, c0, c1), c2);
}

ckks::Ctxt ckks::rescale(const Ctxt & c, Integer idelta)
{
    int lv = c.level - 1;
    if (lv < 0) nevers("negative level");
    Poly c0 = rescaleRound(c.c0, idelta);
    Poly c1 = rescaleRound(c.c1, idelta);
    return Ctxt(lv, c0, c1);
}

ckks::Ctxt ckks::mulExt(const Ctxt & a, const Ctxt & b, const Param & p, const EkExt & ek)
{
    Ctxt3 c3 = mul3(a, b, p);
    Ctxt c2 = relinExt(c3, p, ek);
    Ctxt c2sc = rescale(c2, p.penc.idelta);
    return c2sc;
}

ckks::Ctxt ckks::relinExt(const Ctxt3 & c, const Param & par, const EkExt & ek)
{
    Ctxt r = c; // slice object

    Integer q = par.q(c.level);
    Integer pq = ek.P * q;

    auto d2 = scaleUp(c.c2, q, ek.P, pq);
    auto d2eka = mul(d2, ek.a, pq);
    auto d2ekb = mul(d2, ek.b, pq);
    auto pa = div(d2eka, ek.P, pq);
    auto pb = div(d2ekb, ek.P, pq);
    r.c0 = add(r.c0, pb, q);
    r.c1 = add(r.c1, pa, q); // FIXME pb???
    return r;
}

poly::Poly ckks::genPolyRq(int n, RndStream & rs, Integer q)
{
    Poly r;
    for (int i = 0; i < n; i++) r += rs.getRq(q);
    return r;
}

poly::Poly ckks::genPolyEr(int n, RndStream & rs)
{
    Poly r;
    for (int i = 0; i < n; i++) r += rs.getEr();
    return r;
}

poly::Poly ckks::genPolyR2(int n, RndStream & rs)
{
    Poly r;
    for (int i = 0; i < n; i++) r += Integer(rs.getR2());
    return r;
}

double roundd(double prec, double x)
{
    bool pos = true;
    if (x < 0) { pos = false; x = -x; }
    x /= prec;
    x += 0.5;
    x = (double)(unsigned long long)(x);
    x *= prec;
    return x == 0 ? 0 : pos ? x : -x;
}

std::vector<cx> operator*(const std::vector<cx> & v, double c)
{
    std::vector<cx> r(v);
    for (auto & x : r) x *= c;
    return r;
}

void ckks::ParamEncode::init()
{
    using namespace std::complex_literals;
    int M = n * 2;
    vxi.resize(M);
    static const double pi = std::acos(-1);
    for (int i = 0; i < M; i++)
    {
        //vxi[i] = std::exp(-1i * pi * (1.0 * i) / 4.0);
        vxi[i] = std::exp(1i * pi * (1.0 * i) / 4.0);
    }
}

ckks::Sk::Sk(RndStream & rs, int sz)
{
    n = sz;
    s = genPolyR2(n, rs);
}


Integer ckks::RndStream::getRq(Integer q)
{
    if (1) return Integer(100);
    Integer & a = rq;
    Integer b = ++a;
    b += (b + q / 100) * (q / 100);
    return (b % q);
    ///if (a <= q) a = 0;
    ///return a;
}

int ckks::RndStream::getR2()
{
    if (1) return 1; // FIXME
    int & a = r2;
    return ((++a) % 3 - 1);
}

Integer ckks::RndStream::getEr()
{
    if (1) return Integer(0);
    int & a = er;
    return Integer((++a) % 5 - 2);
}

ckks::Pk::Pk(Sk sk, Param p, RndStream & rs) : n(sk.n)
{
    Poly m0(sk.n, Integer(0));
    Ctxt c = encrypt(sk, m0, p, rs);
    p0 = c.c0;
    p1 = c.c1;
}

void ckks::ParamQx::forceNttValues()
{
    if (ntt::disabled) return;

    for ( auto & x : vql) x = ntt::findCloseQ(penc.n, x);
}

ckks::EkExt::EkExt(Sk sk, Param p, RndStream & rs)
{
    auto Q = p.qL();
    P = Q; // about right - change if need
    P = 2; //* FIXME
    auto PQ = P * p.qL();
    const auto & q = PQ;

    Poly s = sk.s;
    a = genPolyRq(sk.n, rs, q);
    Poly e = genPolyEr(sk.n, rs);
    if (0) e = Poly(sk.n, Integer(0));

    s = rangeUp(s, q);
    e = rangeUp(e, q);

    // b = -a*SK + e + P*SK*SK
    // a = a
    auto x1 = neg(a, q);
    auto x2 = mul(x1, s, q);
    auto x3 = add(x2, e, q);
    auto x4 = mul(s, s, q);
    Poly x5 = mul(x4, P, q);
    b = add(x3, x5, q);
}

string ckks::Param::print() const
{
    std::ostringstream os;
    os << "n = " << penc.n << " \t| polynomial size\n";
    os << "idelta = " << penc.idelta << " \t| Integer delta\n";
    os << "ddelta = " << penc.ddelta << " \t| double delta\n";
    os << "vxi size = " << penc.vxi.size() << " \t| powers of complex root\n";
    os << "levels = " << levels << " \t| Q levels\n";
    os << "values = {";
    for (auto x : vql) os << ' ' << x;
    os << " }" << " \t| Q values\n";
    os << "w = " << w << " \t| Digit decomposition value\n";
    return os.str();
}
