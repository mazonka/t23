#include <iostream> // debug
using std::cout;

#include "ntt.h"
#include "ckkselem.h"

using poly::Poly;
using poly::PolyRns;

const bool D = true;

poly::Poly ckks::encodeP(const Param & p, const std::vector<cx> & v)
{
    std::vector<Integer> m;
    auto vl = v;
    int sz = (int)v.size();
    auto n = p.penc.n;
    if (2 * sz != n) 
        nevers("vector size does not match size in param");

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

    //cout << "AAA1P " << roundv(1e-3,m1) << '\n';

    auto m2 = m1 * (1.0 / n);

    //cout << "AAA2P " << roundv(1e-3, m2) << '\n';

    std::vector<double> m3;
    for (auto x : m2) m3.push_back(x.real());

    //cout << "AAA3P " << roundv(1e-3, m3) << '\n';

    std::vector<double> m4;
    for (auto x : m3) m4.push_back(x * p.penc.ddelta);

    //cout << "AAA4P " << roundv(1e-3, m4) << '\n';

    Poly m5;
    for (auto x : m4)
    {
        m5 += Integer((signed long long)roundd(1.0, x));
        ///cout << " " << x << " ";
    }

    //cout << "AAA5P " << m5 << '\n';

    return { m5 };
}

poly::PolyRns ckks::encodeR(const Param & p, const std::vector<cx> & v, const rns_ns::Rns & rns)
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

    //cout << "AAA1R " << roundv(1e-3, m1) << '\n';

    auto m2 = m1 * (1.0 / n);

    //cout << "AAA2R " << roundv(1e-3, m2) << '\n';

    std::vector<double> m3;
    for (auto x : m2) m3.push_back(x.real());

    //cout << "AAA3R " << roundv(1e-3, m3) << '\n';

    std::vector<double> m4;
    for (auto x : m3) m4.push_back(x * p.penc.ddelta);

    //cout << "AAA4R " << roundv(1e-3, m4) << '\n';

    //std::vector<Integer> m5;
    PolyRns m5(rns);
    for (auto x : m4)
    {
        double val = roundd(1.0, x);
        rns_ns::RnsForm rf(rns, 0, val);
        m5 += rf;
        ///cout << " " << x<< ":" <<val << " "; // AAA
    }

    //cout << "AAA5R " << m5 << '\n';

    return { m5 };
}

std::vector<cx> ckks::decodeP(const Param & p, const Poly & m)
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

std::vector<cx> ckks::decodeR(const Param & p, const PolyRns & m, const rns_ns::Rns & rns)
{
    auto n = p.penc.n;

    ///double dynr = m.rns_ptr->dynrangeDbl();
    ///double dyn2 = dynr / 2;
    std::vector<double> m1;
    ///for (auto x : m.v) m1.push_back(todbl(x));
    {
        int sz = m.polysize();
        for (int i = 0; i < sz; i++)
        {
            rns_ns::RnsForm rf = m.rnsForm(i);
            double x = rf.blendDbl(true);
            ///double x = rf.blendDbl(false);
            ///if (1) // FIXME test
            ///{
            ///    auto y = rf.blend_();
            ///    if (double(y) != x) never;
            ///}
            ///if (x > dyn2) x -= dynr;
            m1.push_back(x);
        }
    }

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

ckks::CtxtP ckks::encryptP(SkP sk, Poly m, Param p, RndStream & rs)
{
    auto q = p.qL_();

    Poly s = sk.s;
    Poly a = genPolyRqP(sk.n, rs, q);
    Poly e = genPolyErP(sk.n, rs);

    s = rangeUpP(s, q);
    e = rangeUpP(e, q);
    auto mR = rangeUpP(m, q);

    // c0 = -a*SK + e + m
    // c1 = a
    auto x1 = neg(a, q);
    auto x2 = mul(x1, s, q);
    auto x3 = add(x2, e, q);
    auto x4 = add(x3, mR, q);

    if(D) cout << "AAA P s=" << s << " a=" << a << " e=" << e << " m=" << mR << '\n';
    if(D) cout << "AAA P x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << '\n';

    CtxtP r(p.levels, x4, a);
    return r;
}

ckks::CtxtR ckks::encryptR(SkR sk, PolyRns mr, Param p, RndStream & rs)
{
    using namespace rns_ns;
    const auto & rns = *mr.rns_ptr;
    RnsForm q = p.qLrns(rns);

    PolyRns s = sk.s;
    PolyRns a = genPolyRqR(sk.n, rs, q, rns);
    PolyRns e = genPolyErR(sk.n, rs, rns);

    if (1) // FIXME16
    {
        e.assign(0, 1);
        e.assign(1, 0);
        a.assign(0, 1000000);
        a.assign(1, 0);
    }

    // no need to range since rns
    //s = rangeUpP(s, q);
    //e = rangeUpP(e, q);
    //auto mR = rangeUpP(m, q);

    // c0 = -a*SK + e + m
    // c1 = a
    auto x1 = neg(a);
    auto x2 = mul(x1, s);
    auto x3 = add(x2, e);
    auto x4 = add(x3, mr);

    if(D) cout << "AAA R s=" << s << " a=" << a << " e=" << e << " m=" << mr << '\n';
    if(D) cout << "AAA R x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << '\n';

    CtxtR r(p.levels, x4, a);
    return r;
}


ckks::CtxtP ckks::encryptP(PkP pk, Poly m, Param p, RndStream & rs)
{
    // c0 = PK0*u+e1+M
    // c1 = PK1*u+e2
    auto q = p.qL_();
    Poly e0 = genPolyErP(pk.n, rs);
    Poly e1 = genPolyErP(pk.n, rs);
    Poly u = genPolyR2P(pk.n, rs);

    u = rangeUpP(u, q);
    e0 = rangeUpP(e0, q);
    e1 = rangeUpP(e1, q);
    auto mR = rangeUpP(m, q);

    auto x1 = mul(pk.p0, u, q);
    auto x2 = add(x1, e0, q);
    auto x3 = add(x2, mR, q);
    auto x4 = mul(pk.p1, u, q);
    auto x5 = add(x4, e1, q);

    //cout << "AAA P e0=" << e0 << " a=" << a << " e=" << e << " m=" << mr << '\n';
    //cout << "AAA P x1=" << x1 << " x4=" << x4 << '\n';

    CtxtP r(p.levels, x3, x5);
    return r;
}

ckks::CtxtR ckks::encryptR(PkR pk, poly::PolyRns mr, Param p, RndStream & rs)
{
    // c0 = PK0*u+e1+M
    // c1 = PK1*u+e2
    //auto q = p.qL_();
    const rns_ns::Rns & rns = *mr.rns_ptr;
    PolyRns e0 = genPolyErR(pk.n, rs, rns);
    PolyRns e1 = genPolyErR(pk.n, rs, rns);
    PolyRns u = genPolyR2R(pk.n, rs, rns);

    //u = rangeUpP(u, q);
    //e0 = rangeUpP(e0, q);
    //e1 = rangeUpP(e1, q);
    //auto mR = rangeUpP(m, q);

    auto x1 = mul(pk.p0, u);
    auto x2 = add(x1, e0);
    auto x3 = add(x2, mr);
    auto x4 = mul(pk.p1, u);
    auto x5 = add(x4, e1);

    CtxtR r(p.levels, x3, x5);
    return r;
}

poly::Poly ckks::decryptP(SkP sk, CtxtP c, Param p)
{
    auto q = p.q_(c.level);

    Poly s = sk.s;

    s = rangeUpP(s, q);
    auto c0 = rangeDownP(c.c0, q);
    auto c1 = rangeDownP(c.c1, q);

    // m = c0 + c1*SK mod q0
    auto x1 = mul(c1, s, q);
    auto x2 = add(c0, x1, q);

    auto m = rangeCenterP(x2, q);

    return m;
}

poly::Poly ckks::decryptP3(SkP sk, Ctxt3P c, Param p)
{
    auto q = p.q_(c.level);

    Poly s = sk.s;

    s = rangeUpP(s, q);
    auto c0 = rangeDownP(c.c0, q);
    auto c1 = rangeDownP(c.c1, q);
    auto c2 = rangeDownP(c.c2, q);

    // m = c0 + c1*SK + c2*SK^2 mod q0
    auto x1 = mul(c1, s, q);
    auto x2 = add(c0, x1, q);
    auto x3 = mul(c2, s, q);
    auto x4 = mul(x3, s, q);
    auto x5 = add(x2, x4, q);

    auto m = rangeCenterP(x5, q);

    return m;
}

poly::PolyRns ckks::decryptR(SkR sk, CtxtR c, Param p)
{
    //auto q = p.q_(c.level);

    PolyRns s = sk.s;

    //s = rangeUpP(s, q);
    //auto c0 = rangeDownP(c.c0, q);
    //auto c1 = rangeDownP(c.c1, q);

    // m = c0 + c1*SK mod q0
    auto x1 = mulDrop(c.c1, s);
    auto x2 = add(c.c0, x1);

    //auto m = rangeCenterP(x2, q);

    //return m;
    return x2;
}

poly::PolyRns ckks::decryptR3(SkR sk, Ctxt3R c, Param p)
{
    PolyRns s = sk.s;

    // m = c0 + c1*SK + c2*SK^2 mod q0
    auto x1 = mul(c.c1, s);
    auto x2 = add(c.c0, x1);
    auto x3 = mul(c.c2, s);
    auto x4 = mul(x3, s);
    auto x5 = add(x2, x4);

    return x5;
}

ckks::CtxtP ckks::add(const CtxtP & a, const CtxtP & b, Param p)
{
    int lv = a.level;
    if (b.level != lv) never;
    auto q = p.q_(lv);
    Poly c0 = poly::add(a.c0, b.c0, q);
    Poly c1 = poly::add(a.c1, b.c1, q);
    return CtxtP(lv, c0, c1);
}

ckks::CtxtR ckks::add(const CtxtR & a, const CtxtR & b)
{
    int lv = a.level;
    if (b.level != lv) never;
    ///auto q = p.q_(lv);
    PolyRns c0 = poly::add(a.c0, b.c0);
    PolyRns c1 = poly::add(a.c1, b.c1);
    return CtxtR(lv, c0, c1);
}

ckks::Ctxt3P ckks::mul3(const CtxtP & a, const CtxtP & b, Param p)
{
    int lv = a.level;
    if (b.level != lv) never;
    auto q = p.q_(lv);

    Poly c0 = mul(a.c0, b.c0, q);
    Poly c1x1 = mul(a.c0, b.c1, q);
    Poly c1x2 = mul(a.c1, b.c0, q);
    Poly c2 = mul(a.c1, b.c1, q);
    Poly c1 = add(c1x1, c1x2, q);
    CtxtP c01(lv, c0, c1);
    return Ctxt3P(c01, c2);
}

ckks::CtxtP ckks::automorphP(const CtxtP& a, Param p)
{
    int lv = a.level;
    auto q = p.q_(lv);

    auto r = a;
    r.c0 = automorph(r.c0);
    r.c1 = automorph(r.c1);
    return r;
}

ckks::CtxtR ckks::automorphR(const CtxtR& a, Param p)
{
    int lv = a.level;
    auto q = p.q_(lv);

    auto r = a;
    r.c0 = automorph(r.c0);
    r.c1 = automorph(r.c1);
    return r;
}

ckks::Ctxt3R ckks::mul3(const CtxtR & a, const CtxtR & b)
{
    int lv = a.level;
    if (b.level != lv) never;

    PolyRns c0 = mul(a.c0, b.c0);
    PolyRns c1x1 = mul(a.c0, b.c1);
    PolyRns c1x2 = mul(a.c1, b.c0);
    PolyRns c2 = mul(a.c1, b.c1);
    PolyRns c1 = add(c1x1, c1x2);
    CtxtR c01(lv, c0, c1);
    return Ctxt3R(c01, c2);
}

static poly::Poly rescaleRoundRecenter(const poly::Poly & pc,
                                       Integer idelta, Integer q)
{
    auto a = poly::rangeCenterP(pc, q);
    Poly b = rescaleRound(a, idelta);
    Poly d = poly::rangeUpP(b, q);
    return d;
}

ckks::Ctxt3P ckks::rescale(const Ctxt3P & c, Integer idelta, Param par)
{
    int lv = c.level - 1;
    if (lv < 0) nevers("negative level");
    Integer q = par.q_(c.level);
    if (0) // old calc
    {
        Poly c0 = rescaleRound(c.c0, idelta);
        Poly c1 = rescaleRound(c.c1, idelta);
        Poly c2 = rescaleRound(c.c2, idelta);
        return Ctxt3P(CtxtP(lv, c0, c1), c2);
    }
    Poly c0 = rescaleRoundRecenter(c.c0, idelta, q);
    Poly c1 = rescaleRoundRecenter(c.c1, idelta, q);
    Poly c2 = rescaleRoundRecenter(c.c2, idelta, q);
    return Ctxt3P(CtxtP(c.level, c0, c1), c2);
}

ckks::Ctxt3R ckks::rescale(const Ctxt3R & c, Integer idelta)
{
    int lv = c.level - 1;
    if (lv < 0) nevers("negative level");
    PolyRns c0 = rescaleRoundRns(c.c0, idelta);
    PolyRns c1 = rescaleRoundRns(c.c1, idelta);
    PolyRns c2 = rescaleRoundRns(c.c2, idelta);
    return Ctxt3R(CtxtR(lv, c0, c1), c2);
}

ckks::CtxtP ckks::rescale(const CtxtP & c, Integer idelta, Param par)
{
    if (0) // old code
    {
        int lv = c.level - 1;
        if (lv < 0) nevers("negative level");
        Poly c0 = rescaleRound(c.c0, idelta);
        Poly c1 = rescaleRound(c.c1, idelta);
        return CtxtP(lv, c0, c1);
    }

    Integer q = par.q_(c.level);
    Poly c0 = rescaleRoundRecenter(c.c0, idelta, q);
    Poly c1 = rescaleRoundRecenter(c.c1, idelta, q);
    return CtxtP(c.level, c0, c1);
}

ckks::CtxtP ckks::rescaleLevelP(const CtxtP & c, Param par)
{
    int level = c.level;
    if (level < 1) nevers("negative level");
    Integer q = par.q_(level);
    Integer idelta = par.vqs[level];
    Poly c0 = rescaleRoundLevel(c.c0, idelta);
    Poly c1 = rescaleRoundLevel(c.c1, idelta);
    return CtxtP(level - 1, c0, c1);
}

ckks::CtxtR ckks::rescale(const CtxtR & c, Integer idelta, Param par)
{
    int lv = c.level - 1;
    if (lv < 0) nevers("negative level");
    PolyRns c0 = rescaleRoundRns(c.c0, idelta);
    PolyRns c1 = rescaleRoundRns(c.c1, idelta);
    return CtxtR(lv, c0, c1);
}

ckks::CtxtR ckks::rescaleLevelR(const CtxtR & c, const rns_ns::RnsShrinkRound & dat)
{
    int lv = c.level;
    if (lv < 1) nevers("negative level");
    ///Integer idelta = par.vqs[lv];
    ///PolyRns c0 = rescaleRoundRnsLevel(c.c0, idelta);
    ///PolyRns c1 = rescaleRoundRnsLevel(c.c1, idelta);
    PolyRns c0 = c.c0.shrink(dat);
    PolyRns c1 = c.c1.shrink(dat);
    return CtxtR(lv - 1, c0, c1);
}

ckks::CtxtP ckks::relinExtP(const Ctxt3P & c, const Param & par, const EkExtP & ek)
{
    CtxtP r = c.slice(); // slice object

    Integer q = par.q_(c.level);
    Integer pq = ek.P * q;

    //auto d2 = scaleUp(c.c2, q, ek.P, pq);
    auto d2 = rangeUpP(c.c2, q);

    auto d2eka = mul(d2, ek.a, pq);
    auto d2ekb = mul(d2, ek.b, pq);
    auto pa = div(d2eka, ek.P, pq);
    auto pb = div(d2ekb, ek.P, pq);
    r.c0 = add(r.c0, pb, q);
    r.c1 = add(r.c1, pa, q);
    if(D) cout << "AAA " << __func__ << " d2=" << d2 << '\n'
         << " ek:a:b=" << ek.a << ek.b << '\n'
         << " d2ek=" << d2eka << d2ekb << '\n';
    return r;
}

ckks::CtxtR ckks::relinExtR(const Ctxt3R & c, const Param & par, const EkExtR & ek)
{
    CtxtR r = c.slice(); // slice object

    ///never;
    ///Integer q = par.q_(c.level);
    ///Integer pq = ek.P * q;

    //auto d2 = scaleUp(c.c2, q, ek.P, pq);
    const auto & rext = *ek.a.rns_ptr;
    auto d2 = c.c2.rebase(rext);
    auto d2eka = mul(d2, ek.a);
    auto d2ekb = mul(d2, ek.b);
    ///auto pa = div(d2eka, ek.P, pq);
    ///auto pb = div(d2ekb, ek.P, pq);
    auto pa = d2eka.shrink(ek.rshrink);
    auto pb = d2ekb.shrink(ek.rshrink);
    r.c0 = add(r.c0, pb);
    r.c1 = add(r.c1, pa);
    if(D) cout << "AAA " << __func__ << " d2=" << d2 << '\n'
         << " ek:a:b=" << ek.a << ek.b << '\n'
         << " d2ek=" << d2eka << d2ekb << '\n';
    return r;
}

ckks::CtxtP ckks::mulExtP(const CtxtP & a, const CtxtP & b, const Param & p, const EkExtP & ek)
{
    Ctxt3P c3 = mul3(a, b, p);
    CtxtP c2 = relinExtP(c3, p, ek);
    CtxtP c2sc = rescaleLevelP(c2, p);
    return c2sc;
}

ckks::CtxtR ckks::mulExtR(const CtxtR & a, const CtxtR & b, const Param & p, const EkExtR & ek, const rns_ns::RnsShrinkRound & datQ)
{
    Ctxt3R c3 = mul3(a, b);
    CtxtR c2 = relinExtR(c3, p, ek);
    CtxtR c2sc = rescaleLevelR(c2, datQ);
    return c2sc;
}


poly::Poly ckks::genPolyRqP(int n, RndStream & rs, Integer q)
{
    Poly r;
    for (int i = 0; i < n; i++) r += getRq(rs, q);
    return r;
}

poly::Poly ckks::genPolyErP(int n, RndStream & rs)
{
    Poly r;
    for (int i = 0; i < n; i++) r += rs.getEr();
    return r;
}

poly::Poly ckks::genPolyR2P(int n, RndStream & rs)
{
    Poly r;
    for (int i = 0; i < n; i++) r += Integer(rs.getR2());
    return r;
}

poly::PolyRns ckks::genPolyRqR(int n, RndStream & rs, rns_ns::RnsForm q, const rns_ns::Rns & rns)
{
    PolyRns r(rns);
    ///for (int i = 0; i < n; i++) r += rs.getRq(q);
    for (int i = 0; i < n; i++) r += getRqRns(rs, q);
    return r;
}

poly::PolyRns ckks::genPolyErR(int n, RndStream & rs, const rns_ns::Rns & rns)
{
    PolyRns r(rns);
    for (int i = 0; i < n; i++)
        r += rns_ns::RnsForm(rns, rs.getEr());
    return r;
}

poly::PolyRns ckks::genPolyR2R(int n, RndStream & rs, const rns_ns::Rns & rns)
{
    PolyRns r(rns);
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

ckks::SkP::SkP(RndStream & rs, int sz)
{
    n = sz;
    s = genPolyR2P(n, rs);
}

ckks::SkR::SkR(RndStream & rs, int sz, const rns_ns::Rns & rns):
    n(sz), s(genPolyR2R(sz, rs, rns))
{
    ///n = sz;
    ///s = genPolyR2R(n, rs);
    if (1) // FIXME16
    {
        s.assign(0, 0);
        s.assign(1, 0);
    }
}

Integer ckks::RndStream::getRqSeed()
{
    if (0) return Integer(0);
    Integer & a = rq;
    Integer b = ++a;
    return b;
}

Integer ckks::getRq(RndStream & rs, Integer q0)
{
    //if (0) return Integer(0);
    //Integer b = rs.getRqSeed();
    //auto q = q0;

    //if (0) // new version
    //{
    //    q -= 1;
    //    q /= 8;
    //    q += b;
    //    b += b * q;
    //    b = b % q0;
    //    return b;
    //}
    //else // old version
    //{
    //    Integer b = rs.getRqSeed();
    //    b += (b + q / 100) * (q / 100);
    //    return (b % q);
    //}

    if (0) return Integer(0); // FIXME
    Integer b0 = rs.getRqSeed();
    //b += (b + q0 / 100) * (q0 / 10);
    //return (b % q0);

    //if (b0 == 8) return 0;
    //if (b0 == 7) return 0; // 217146815263
    //if (b0 == 6) return 0;
    //if (b0 == 5) return 0;
    //if (b0 == 4) return 0;
    //if (b0 == 3) return 0;
    //if (b0 == 2) return 0+0*263678;
    //if (b0 == 1) return 1+0*131838;
    //return q0-1;


    if (1) // new version
    {
        auto q1 = q0;
        auto b1 = b0;
        q1 -= 1;
        q1 /= 8;
        q1 += b1;
        b1 += b1 * q1;
        if (b1 < 0) never;
        b1 = b1 % q0;
        b1 = q0 - b1; // addition 01
        std::ignore = b1;
        return b1;
    }
    else // old version
    {
        auto q2 = q0;
        auto b2 = b0;
        b2 += (((b2 + q2 / 100) % q0) * (q2 / 10)) % q0;
        b2 %= q0;
        b2 += b2 / 1324; // 1324-bad 1325-ok - poly 4
        b2 %= q0;
        //b2 += (b2 * 2000)%q0; // poly 2
        if (b2 < 0) never;
        b2 %= q0;
        return b2;
    }

}

rns_ns::RnsForm ckks::getRqRns(RndStream & rs, const rns_ns::RnsForm & fqm)
{
    using rns_ns::RnsForm;

    if (0) // old code
    {
        // FIXME this function is a mess
        // need simple algorithm

        //auto vs = fq.values();
        //Integer q0(0);
        //if (fq.islowval() && fq.lowval() == 0)
        //    for (int i = 0; i < fq.size(); i++)
        //        q0 = vs[i] = fq.rns()->q(i);

        //vint vr;
        //if (q0 == 0)
        //{
        //    Integer sum = 0;
        //    for (auto q : vs) sum += q;
        //    auto b = rs.getRq(sum);
        //    for (auto q : vs) vr.push_back(b);
        //    return rns_ns::RnsForm(fq.rns(), vr);
        //}

        //// this is incorrect, but we need only one call to getRq
        //// to match non-rns computation. It can be changed later.
        //return rns_ns::RnsForm(fq.rns(), rs.getRq(q0));
    }

    auto rns = fqm.rns();
    if (0) return RnsForm(rns, 0); // FIXME
    Integer b = rs.getRqSeed();
    RnsForm fb(rns, b), fq { fqm };
    //fq.setM();
    //fb += (fb + fq / Integer{ 100 }) * (fq / Integer{ 10 });
    //return (fb % fq);

    fq -= Integer { 1 };
    ///fq.blend_();
    fq /= 8;
    ///fq.blend_();
    fq += fb;
    ///fq.blend_();
    fb += fb * fq;
    ///fb.blend_();
    //fb = fb % q0;
    if (0) ///
    {
        auto xb = fb.blend_();
        auto xd = rns->dynrange_();
        std::ignore = xb;
        std::ignore = xd;
    }
    fb = - fb; // addition 01

    return fb;
}

int ckks::RndStream::getR2()
{
    if (0) return 0; // FIXME
    int & a = r2;
    ++a;

    if (a == 1) return 1;
    return -1;

    auto r = ((a + 1) % 3 - 1);
    return r;
}

Integer ckks::RndStream::getEr()
{
    if (0) return Integer(0);
    //return Integer(0); // FIXME *************************************
    int & a = er;
    ++a;

    //if (a == 0) return Integer(0); // FIXME *************************************
    //if (a == 1) return Integer(0); // FIXME *************************************
    //if (a == 2) return Integer(0); // FIXME *************************************
    //if (a == 3) return Integer(1); // FIXME *************************************
    //if (a == 4) return Integer(0); // FIXME *************************************
    //if (a == 5) return Integer(0); // FIXME *************************************
    //if (a == 6) return Integer(0); // FIXME *************************************
    //if (a == 7) return Integer(0); // FIXME *************************************
    //if (a == 8) return Integer(0); // FIXME *************************************
    //if (a == 9) return Integer(0); // FIXME *************************************
    //if (a == 10) return Integer(0); // FIXME *************************************
    //if (a == 11) return Integer(0); // FIXME *************************************
    return Integer(a % 5 - 2);
}

ckks::PkP::PkP(SkP sk, Param p, RndStream & rs)
    : n(sk.n)
{
    Poly m0(sk.n, Integer(0));
    CtxtP c = encryptP(sk, m0, p, rs);
    p0 = c.c0;
    p1 = c.c1;
}

ckks::PkR::PkR(SkR sk, Param p, RndStream & rs)
    : n(sk.n), p0(sk.s.rns_ptr), p1(sk.s.rns_ptr)
{
    PolyRns m0(sk.s.rns_ptr, sk.n, Integer(0));
    CtxtR c = encryptR(sk, m0, p, rs);
    p0 = c.c0;
    p1 = c.c1;
}

Integer ckks::ParamQx::qL_() const
{
    Integer r(1);
    for (auto x : vqs) r *= x;
    return r;
}

Integer ckks::ParamQx::q_(int l) const
{
    Integer r(1);
    for (int i = 0; i < l + 1; i++) r *= vqs[i];
    return r;
}

rns_ns::RnsForm ckks::ParamQx::qLrns(const rns_ns::Rns & rns) const
{
    using rns_ns::RnsForm;
    RnsForm r(rns, 1);
    for (auto q : vqs) r *= RnsForm(rns, q);
    return r;
}

std::vector<Integer> ckks::ParamQx::findNttCoprimes(Integer q)
{
    return ntt::findCoprimes(penc.n, q, levels + 1);
}

ckks::EkExtP::EkExtP(SkP sk, Param p, RndStream & rs, Integer newp)
{
    auto Q = p.qL_();

    if (newp == 0)
        P = Q; // about right - change if need
    else
        P = newp;

    auto PQ = P * Q;
    const auto & q = PQ;

    Poly s = sk.s;
    a = genPolyRqP(sk.n, rs, q);
    Poly e = genPolyErP(sk.n, rs);
    if (0) e = Poly(sk.n, Integer(0));

    s = rangeUpP(s, q);
    e = rangeUpP(e, q);

    // b = -a*SK + e + P*SK*SK
    // a = a
    auto x1 = neg(a, q);
    auto x2 = mul(x1, s, q);
    auto x3 = add(x2, e, q);
    auto x4 = mul(s, s, q);
    Poly x5 = mul(x4, P, q);
    b = add(x3, x5, q);

    if(D) cout << "AAA " << __func__ << " PQ=" << q << " s=" << s << " a=" << a << " e=" << e << " b=" << b << '\n';
    if(D) cout << "AAA2 " << "x1x2x3x4x5 " << x1 << x2 << x3 << x4 << x5 << '\n';
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
    for (auto x : vqs) os << ' ' << x;
    os << " }" << " \t| Q values\n";
    os << "w = " << w << " \t| Digit decomposition value\n";
    return os.str();
}

ckks::EkExtR::EkExtR(SkR sk, Param p, RndStream & rs, rns_ns::Rns & rext,
                     rns_ns::RnsShrinkRound rshr)
    : P { 0 }, b(sk.s.rns_ptr), a(sk.s.rns_ptr), rshrink(rshr)
{
    //auto Q = p.qL_();
    //P = Q; // about right - change if need
    //auto PQ = P * p.qL_();
    //const auto& q = PQ;

    //rns_ns::RnsForm qrf(rext, 0, p.vqs);
    rns_ns::RnsForm qrf(rext, 0); // all range

    PolyRns s = sk.s;
    a = genPolyRqR(sk.n, rs, qrf, rext);
    PolyRns e = genPolyErR(sk.n, rs, rext);

    ///s = rangeUpP(s, q);
    ///e = rangeUpP(e, q);


    //auto se = s.rebase(rext);
    auto se = s.modswap(rext);
    auto PinPQ = rshrink.PinQ.rebaseAdd(rext);

    // b = -a*SK + e + P*SK*SK
    // a = a
    auto x1 = neg(a);
    auto x2 = mul(x1, se);
    auto x3 = add(x2, e);
    auto x4 = mul(se, se);
    PolyRns x5 = mul(x4, PinPQ);
    b = add(x3, x5);

    if(D) cout << "AAA " << __func__ << " PQ=" << rext.dynrange_() << " s=" << se << " a=" << a << " e=" << e << " b=" << b << '\n';
    if(D) cout << "AAA2 " << "x1x2x3x4x5 " << x1 << x2 << x3 << x4 << x5 << '\n';
}
