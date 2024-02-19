///#include <iostream> /// debug
///using std::cout;

#include <set>

#include "integer.h"
#include "poly.h"
#include "ntt.h"

///rns_ns::RnsMrs poly::PolyRns::rns;

std::ostream & operator<<(std::ostream & os, const poly::Poly & p)
{
    os << "{";
    for (auto x : p.v) os << ' ' << x;
    os << " }";
    return os;
}


poly::Poly poly::rangeUpP(Poly a, Integer q)
{
    Poly r = a;
    for (auto & x : r.v) while ( x < 0 ) x += q;
    return r;
}

poly::Poly poly::rangeDownP(Poly a, Integer q)
{
    Poly r = a;
    for (auto & x : r.v) x = x % q;
    return r;
}

poly::Poly poly::rangeCenterP(Poly a, Integer q)
{
    Poly r = a;
    auto q2 = q / 2;
    for (auto & x : r.v) if ( x > q2 ) x -= q;
    return r;
}

poly::Poly poly::neg(Poly a, Integer q)
{
    Poly r;
    for (auto x : a.v) r += q - x;
    return r;
}

poly::PolyRns poly::neg(PolyRns a)
{
    PolyRns r(a);
    r.negInplace();
    return r;
}

poly::PolyRns poly::add(PolyRns a, PolyRns b)
{
    int sz = a.ntows();
    if (sz != b.ntows()) never;

    for (int i = 0; i < sz; i++)
    {
        auto & ati = a.towers[i];
        ati = add(ati, b.towers[i], a.rns_ptr->q(i));
    }
    return a;
}

poly::Poly poly::mulmod_simple(const Poly & a, const Poly & b, Integer q)
{
    int cntr = 0;

    auto n = b.size();  if (a.size() != n) never;
    Poly r(n * 2, Integer(0));

    auto assertRange = [&](Integer x)
    {
        //cout << (++cntr) << " " << x << '\n';
        if (x < Integer(0)) never;
        if (x >= q) never;
    };

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            auto & rij = r.v[i + j];
            ///cout << (++cntr) << " Rij a b : " << rij <<' '<< a.v[i]<<' '<< b.v[j] << '\n';
            rij += modmul(a.v[i], b.v[j], q);
            ///cout << (++cntr) << " Rij : " << rij << '\n';
            if (rij >= q) rij -= q;
            //assertRange(rij);
        }

    for (int i = 0; i < n; i++)
    {
        auto & d = r.v[i];
        auto & u = r.v[i + n];
        if (d < u) d += q;
        d -= u;
        //assertRange(d);
    }
    r.v.resize(n);

    ///cout << "c " << r << '\n';
    return r;
}

poly::Poly poly::mul_simple(const Poly & a, const Poly & b)
{
    auto n = b.size();  if (a.size() != n) never;
    Poly r(n * 2, Integer(0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            auto & rij = r.v[i + j];
            rij += a.v[i] * b.v[j];
        }

    for (int i = 0; i < n; i++)
    {
        auto & d = r.v[i];
        auto & u = r.v[i + n];
        d -= u;
    }
    r.v.resize(n);

    return r;
}

poly::Poly poly::mulmod_elwise(const Poly & a, const Poly & b, Integer q)
{
    auto sz = b.size();  if (a.size() != sz) never;
    Poly r(sz);
    for (int i = 0; i < sz; i++) r.v[i] = modmul(a.v[i], b.v[i], q);
    return r;
}

poly::Poly poly::mulmod_btrfly(const Poly & a, const Poly & b, Integer q)
{
    if (ntt::disabled) return mulmod_simple(a, b, q);

    using namespace ntt;
    auto an = nttBfly(a, q);
    auto bn = nttBfly(b, q);
    auto cn = mulmod_elwise(an, bn, q);
    auto c = ittBfly(cn, q);
    return c;
}


poly::Poly poly::add(Poly a, Poly b, Integer q)
{
    auto sz = b.size();  if (a.size() != sz) never;

    Poly r(sz);
    for (int i = 0; i < sz; i++)
    {
        r.v[i] = a.v[i] + b.v[i];
        if (r.v[i] >= q) r.v[i] -= q;
        if (1) if (r.v[i] >= q) nevers("input overflow");
    }
    return r;
}

// requires recenterd input and outputs recentered
poly::Poly poly::rescaleRound(const Poly & a, Integer idelta)
{
    auto d2 = idelta / 2;
    Poly r(a);
    for (auto & x : r.v)
    {
        if (x < 0)
            x = -((-x + d2) / idelta);
        else
            x = (x + d2) / idelta;
    }
    return r;
}

poly::Poly poly::rescaleRoundLevel(const Poly & a, Integer idelta)
{
    auto d2 = idelta / 2;
    Poly r(a);
    for (auto & x : r.v)
    {
        if (x < 0) never;
        //    x = -((-x + d2) / idelta);
        //else
        x = (x + d2) / idelta;
    }
    return r;
}


//poly::Poly poly::rescaleFloor(const Poly & a, Integer w)
//{
//    Poly r(a);
//    for (auto & x : r.v) x /= w;
//    return r;
//}

poly::Poly poly::mul(Poly a, Integer b, Integer q)
{
    auto sz = a.size();

    Poly r(sz);
    for (int i = 0; i < sz; i++)
        r.v[i] = modmul(a.v[i], b, q);

    return r;
}

poly::Poly poly::div(Poly a, Integer b, Integer q)
{
    a = rescaleRound(a, b);
    a = rangeDownP(a, q);
    return a;
}

poly::Poly poly::scaleUp(const Poly & a, Integer Q, Integer P, Integer PQ)
{
    Poly r = rangeCenterP(a, Q);
    for (auto & x : r.v) x = x * P;
    r = rangeUpP(r, PQ);
    return r;
}


// Poly RNS

void poly::PolyRns::operator+=(Integer x)
{
    rns_ns::RnsForm rf(rns_ptr, x);
    *this += rf;
}

void poly::PolyRns::operator+=(rns_ns::RnsForm x)
{
    if (!x.match(rns_ptr)) never;
    int intows = rns_ptr->size();
    if (towers.empty())
        towers.resize(intows);
    else if (intows != (int)towers.size()) never;

    for (int i = 0; i < intows; i++)
        towers[i] += x.values()[i];
}

poly::PolyRns::PolyRns(const rns_ns::Rns * r, int n, Integer i)
    : rns_ptr(r)
{
    Tower t(n, i);
    for (int i = 0; i < r->size(); i++) towers.push_back(t);
}

poly::PolyRns::Tower poly::PolyRns::getCollapsed(bool forceBlend) const
{
    Poly r;
    int sz = polysize();
    for (int i = 0; i < sz; i++)
    {
        ///int intows = (int).size();
        vint vf;
        for (int j = 0; j < ntows(); j++)
            vf.push_back(towers[j].v[i]);
        rns_ns::RnsForm rf(rns_ptr, vf);
        if (forceBlend)
            r += rf.blend_();
        else
        {
            if (!rf.islowval()) nevers("big value to fit");
            r += rf.lowval();
        }
    }
    return r;
}

rns_ns::RnsForm poly::PolyRns::rnsForm(int i) const
{
    int nt = ntows();
    vint v(nt);
    for (int j = 0; j < nt; j++)
        v[j] = towers[j].v[i];

    return rns_ns::RnsForm(rns_ptr, v);
}

void poly::PolyRns::negInplace()
{
    for (int i = 0; i < ntows(); i++)
    {
        auto q = rns_ptr->q(i);
        towers[i] = neg(towers[i], q);
    }
}

bool poly::PolyRns::match(const PolyRns & b) const
{
    if (rns_ptr != b.rns_ptr) return false;
    if (ntows() != b.ntows()) return false;
    return true;
}

poly::PolyRns poly::PolyRns::rebase(const rns_ns::Rns & nr) const
{
    PolyRns r(nr);
    int sz = polysize();

    for (int i = 0; i < sz; i++)
        r += rnsForm(i).rebaseAny(nr);

    return r;
}

// same as rebase but treats negatives
poly::PolyRns poly::PolyRns::modswap(const rns_ns::Rns & nr) const
{
    PolyRns r(nr);
    int sz = polysize();

    for (int i = 0; i < sz; i++)
        r += rnsForm(i).baseSwap(nr);

    return r;
}

poly::PolyRns poly::PolyRns::shrink(const rns_ns::RnsShrinkRound & rshr) const
{
    const rns_ns::Rns & nr = rshr.Q;
    PolyRns r(nr);
    int sz = polysize();

    for (int i = 0; i < sz; i++)
        r += rnsForm(i).rebaseShrinkRound(rshr);

    return r;
}

poly::PolyRns poly::mul(const PolyRns & a, const PolyRns & b)
{
    int ntows = a.ntows();
    if (ntows != b.ntows()) never;
    auto rns = a.rns_ptr;
    if (rns != b.rns_ptr) never;

    PolyRns ret(a.rns_ptr);
    for (int i = 0; i < ntows; i++)
        ret.towers.emplace_back(mul(a.towers[i], b.towers[i], rns->q(i)));

    return ret;
}

// same as mul but allows to drop excessive towersw in b
poly::PolyRns poly::mulDrop(const PolyRns & a, const PolyRns & b)
{
    int antows = a.ntows();
    int bntows = b.ntows();
    auto arns = a.rns_ptr;
    auto brns = b.rns_ptr;
    if (antows == bntows)
    {
        if (arns != brns) never;
        return mul(a, b);
    }

    PolyRns bdrop(arns);

    auto avqs = arns->getQs();
    std::set<Integer> savqs;
    savqs.insert(avqs.begin(), avqs.end());
    for (int i = 0; i < bntows; i++)
    {
        auto q = brns->getQs()[i];
        if (savqs.find(q) == savqs.end()) continue;
        bdrop.towers.push_back(b.towers[i]);
    }

    if (antows != bdrop.ntows()) never;

    return mul(a, bdrop);
}

poly::PolyRns poly::mul(const PolyRns & a, rns_ns::RnsForm b)
{
    int ntows = a.ntows();
    PolyRns r(a.rns_ptr);
    vint bv = b.values();
    const rns_ns::Rns * rns = a.rns_ptr;
    if (rns != b.rns() ) never;
    vint qs = rns->getQs();
    for (int i = 0; i < ntows; i++)
    {
        Poly p = mul(a.towers[i], bv[i], qs[i]);
        r.towers.push_back(p);
    }
    return r;
}

//poly::PolyRns poly::mul(PolyRns a, Integer b)
//{
//    never;
//    ///return PolyRns();
//}
//
//poly::PolyRns poly::mul(PolyRns a, const vint & vb)
//{
//    never;
//    ///return PolyRns();
//}

poly::PolyRns poly::mul_simple(const PolyRns & a, const PolyRns & b)
{
    if (!a.match(b)) never;
    PolyRns r(a.rns_ptr);
    for (int i = 0; i < a.ntows(); i++)
        r.towers.emplace_back(mul_simple(a.towers[i], b.towers[i]));
    return r;
}

poly::PolyRns poly::rescaleRoundRns(const PolyRns & a, Integer idelta)
{
    ///never;
    ///return PolyRns();
    using rns_ns::RnsForm;
    auto d2 = idelta / 2;
    PolyRns r(a.rns_ptr);
    int sz = a.polysize();

    for (int i = 0; i < sz; i++)
    {
        RnsForm rf = a.rnsForm(i);
        double x = rf.blendDbl(true);
        if (x < 0)
            x = -((-x + d2) / idelta);
        else
            x = (x + d2) / idelta;
        r += RnsForm(*a.rns_ptr, 0, x);
    }
    return r;
}

