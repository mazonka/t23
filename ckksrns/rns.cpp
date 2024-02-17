#include <iostream> /// debug
using std::cout;

#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>

#include "rns.h"
#include "err.h"


std::pair<Integer, int> rns_ns::Rns::blend_(const vint & x) const
{
    int sz = size();
    auto d = dynrange_();

    Integer sum(0);
    int rank = 0;

    for ( int i = 0; i < sz; i++ )
    {
        auto xi = (x[i] * us[i]) % qs[i];
        auto mi = xi * (d / qs[i]);
        sum += mi;
        if (sum >= d)
        {
            sum -= d;
            ++rank;
        }
        if (sum >= d) never;
    }
    return { sum, rank };
}

double rns_ns::Rns::blendDbl(const vint & x, bool recenter) const
{
    int sz = size();
    ///auto d = dynrange_();

    double sum(0);

    double d = 1;
    std::vector<double> mi(sz, 1);
    for (int i = 0; i < sz; i++)
    {
        d *= qs[i];
        for (int j = 0; j < sz; j++)
            if (i != j) mi[j] *= qs[i];
    }

    for (int i = 0; i < sz; i++)
    {
        double dxi = double(x[i]);
        double dusi = double(us[i]);
        double dqsi = double(qs[i]);
        double xi = std::fmod(dxi * dusi, dqsi);

        sum += mi[i] * xi;
        if (sum >= d) sum -= d;
    }
    if (recenter && sum > d / 2) sum -= d;
    return sum;
}


vint rns_ns::Rns::neg(vint v) const
{
    int sz = size();
    if (::isize(v) != sz) never;
    for (int i = 0; i < sz; i++)
    {
        if (v[i] == 0) continue;
        v[i] = qs[i] - v[i];
    }
    return v;
}

vint rns_ns::Rns::split(Integer x) const
{
    if (x < 0) return neg(split(-x));
    assertRnsInited();
    ///if (x == 0) return vint(qs.size(), 0);
    vint r;
    ///bool bn = false;
    ///if (x < 0)
    ///{
    ///bn = true;
    ///x = -x;
    ///return negate(split(-x));
    ///}

    for (auto q : qs) r.push_back(x % q);
    ///if (bn) negate(r);
    return r;
}

///if (neg)
///r.push_back(q - (x % q));
///else

vint rns_ns::Rns::splitFmod(double x) const
{
    if (x < 0) return neg(splitFmod(-x));
    assertRnsInited();
    vint r;
    for (auto q : qs)
    {
        double a = std::fmod(x, q);
        Integer b = (long long int)(a);
        r.push_back(b);
    }
    return r;
}

Integer rns_ns::Rns::lowval(const vint & v) const
{
    if (!islowval(v)) nevers("rns overflow: cannot represent in q0 range");
    return v[size() - 1];
}

bool rns_ns::Rns::islowval(const vint & v) const
{
    int sz = ::isize(v);
    for (int i = 1; i < sz; i++) // 1
        if (v[0] != v[i]) goto chkmax;
    return true;

chkmax:
    int sz1 = sz - 1;
    for (int i = 0; i < sz1; i++)
        if (v[sz1] % qs[i] != v[i] ) return false;
    return true;
}

void rns_ns::Rns::initL(std::initializer_list<Integer> li)
{
    initV(li);
}

void rns_ns::Rns::initV(const vint & ve)
{
    qs = ve;
    if (qs.empty()) never;
    std::sort(qs.begin(), qs.end());
    // init mu
    for (auto q : qs)
    {
        Integer m = 1, m_ = 1;
        for (auto s : qs)
            if (s != q)
            {
                m = modmul(m, s, q);
                m_ *= s;
            }
        Ms_.push_back(m_);
        us.push_back(modinv(m, q));
    }

    maxPow2digit = -1;
    auto x = qs[0];
    while (x) { ++maxPow2digit; x >>= 1; }
}

string rns_ns::Rns::print() const
{
    std::ostringstream os;
    os << "Dynrange = " << dynrange_() << '\n';
    os << "qs ="; for (auto x : qs) os << ' ' << x; os << '\n';
    os << "Ms ="; for (auto x : Ms_) os << ' ' << x; os << '\n';
    os << "us ="; for (auto x : us) os << ' ' << x; os << '\n';
    os << "maxPow2Digit = " << maxPow2digit << '\n';
    return os.str();
}

Integer rns_ns::Rns::dynrange_() const
{
    Integer r = 1;
    for (auto x : qs) r *= x;
    return r;
}

double rns_ns::Rns::dynrangeDbl() const
{
    double r = 1;
    for (auto x : qs) r *= double(x);
    return r;
}

vint rns_ns::Rns::qsetOp(const vint & a, Op op, const vint & b)
{
    std::set<Integer> s;
    if (op == minus)
    {
        s.insert(a.begin(), a.end());
        for (auto x : b)
        {
            if (s.find(x) == s.end()) nevers("minus must be a subset");
            s.erase(x);
        }
        vint r;
        r.insert(r.end(), s.begin(), s.end());
        return r;
    }

    if (op == plus)
    {
        s.insert(a.begin(), a.end());
        for (auto x : b)
        {
            if (s.find(x) != s.end()) nevers("plus must be exclusive");
            s.insert(x);
        }
        vint r;
        r.insert(r.end(), s.begin(), s.end());
        return r;
    }

    never;
}

rns_ns::RnsForm::RnsForm(const Rns & r, int zero, vint vs) : prns(&r)
{
    RnsForm a(r, 1);
    for (auto x : vs) a *= RnsForm(r, x);
    v = a.values();
}

bool rns_ns::RnsForm::isnegative() const
{
    auto x = *this * Integer { 2 };
    if (x < *this) return true;
    return false;
}

rns_ns::RnsForm rns_ns::RnsForm::invert() const
{
    RnsForm r(*this);
    int sz = size();
    vint qs = prns->getQs();
    for (int i = 0; i < sz; i++)
    {
        if (v[i] == 0) nevers("inversion fails on 0");
        r.v[i] = modinv(v[i], qs[i]);
    }
    return r;
}

rns_ns::RnsForm rns_ns::RnsForm::negate() const
{
    RnsForm r = *this;
    r.v = prns->neg(v);
    return r;
}

// set 0 to M
void rns_ns::RnsForm::setM()
{
    if (!islowval()) return;
    if (lowval() != Integer { 0 }) return;
    v = prns->getQs();
}

int rns_ns::RnsForm::operator+=(const vint & b)
{
    const auto & qs = prns->qs;
    int sz = prns->size();
    int overflows = 0;
    for (int i = 0; i < sz; i++)
    {
        auto [x, y] = modaddOver(v[i], b[i], qs[i]);
        v[i] = x;
        overflows += y;
    }
    return overflows;
}

rns_ns::RnsForm & rns_ns::RnsForm::operator+=(const rns_ns::RnsForm & b)
{
    if (!prns) never;
    const auto & qs = prns->qs;
    int sz = prns->size();
    for (int i = 0; i < sz; i++) v[i] = modadd(v[i], b.v[i], qs[i]);
    return *this;
}

std::pair<rns_ns::RnsForm, rns_ns::RnsForm> rns_ns::RnsForm::divABRQ(const RnsForm & a, const RnsForm & b) const
{
    // slow
    if (a < b)
        return { a, RnsForm(a.prns, 0) };

    if (b.islowval() && b.lowval() == Integer { 0 })
        nevers("division by 0");

    return divABRQ_rec(a, b);
}

std::pair<rns_ns::RnsForm, rns_ns::RnsForm> rns_ns::RnsForm::divABRQ_rec(const RnsForm & a, const RnsForm & b) const
{
    // we assume that b<=a

    if (1) // FIXME turn off
        if (a < b) never;

    RnsForm f0(a.prns, 0), f1(a.prns, 1), f2(a.prns, 2);
    Integer i2 { 2 };

    if ( a == b ) return { f0, f1 };

    auto c = a - b;
    if ( c == b ) return { f0, f2 };
    if ( c < b )  return { c, f1 };

    auto b2 = b * i2;
    auto [R, Q] = divABRQ_rec(a, b2);
    Q *= i2;
    if (R >= b)
    {
        R -= b;
        Q += Integer { 1 };
    }
    return { R, Q };
}

bool rns_ns::RnsForm::operator<(const RnsForm & b) const
{
    vint am = prns->mrs(v);
    vint bm = prns->mrs(b.v);
    int sz = (int)am.size();
    for (int i = sz - 1; i >= 0; i--)
    {
        if (am[i] < bm[i]) return true;
        if (bm[i] < am[i]) return false;
    }
    return false;
}

rns_ns::RnsForm rns_ns::RnsForm::operator-() const
{
    if (!prns) never;
    const auto & qs = prns->qs;
    int sz = prns->size();
    RnsForm r(*this);
    for (int i = 0; i < sz; i++) r.v[i] = qs[i] - v[i];
    return r;
}

rns_ns::RnsForm rns_ns::RnsForm::rebaseAny(const Rns & nr) const
{
    vint mr = prns->mrs(v); // MRS of this value
    vint nqs = prns->getQs();

    RnsForm r(nr, 0);
    for (int i = (int)nqs.size() - 1; i >= 0; i--)
    {
        r *= RnsForm(nr, nqs[i]);
        r += RnsForm(nr, mr[i]);
    }
    return r;
}

rns_ns::RnsForm rns_ns::RnsForm::baseSwap(const Rns & newRns) const
{
    bool neg = isnegative();
    auto r = *this;
    if (neg) r = r.negate();
    r = r.rebaseAny(newRns);
    if (neg) r = r.negate();
    return r;
}

rns_ns::RnsForm rns_ns::RnsForm::rebaseCut(const Rns & nr) const
{
    vint vqs = prns->getQs();
    vint nqs = nr.getQs();
    const vint & xis = v;

    int tsz = prns->size();
    int sz = nr.size();
    vint selectedxis(sz);
    for (int i = 0, j = 0; i < sz; i++, j++)
    {
        while (nqs[i] != vqs[j])
            if (++j == tsz)
                nevers("no match rns");
        selectedxis[i] = xis[j];
    }

    return RnsForm(nr, selectedxis);
}

rns_ns::RnsForm rns_ns::RnsForm::rebaseShrinkRound(const RnsShrinkRound & d) const
{
    /*
    * ((X/P)) = [(X+P/2)/P] = |P^-1_q(X+P/2-|X+P/2|_p)|_q =
    * | P^-1_q ( X_q + |P/2|_q - [ |X_p+P/2|_p ]_q ) |_q
    * constants: P^-1_q=P1q, P/2_q=P2q, P/2=|P/2|_p= P2
    */

    RnsForm X_q = rebaseCut(d.Q);
    RnsForm X_p = rebaseCut(d.P);
    auto a = X_q + d.P2q;
    auto b1 = X_p + d.P2p;
    auto b2 = b1.rebaseAny(d.Q);
    auto c = a - b2;
    auto e = d.P1q * c;

    return e;
}

rns_ns::RnsForm & rns_ns::RnsForm::operator*=(const rns_ns::RnsForm & b)
{
    if (!prns) never;
    const auto & qs = prns->qs;
    int sz = prns->size();
    for (int i = 0; i < sz; i++) v[i] = modmul(v[i], b.v[i], qs[i]);
    return *this;
}

void rns_ns::RnsMrs::myinit()
{
    int sz = size();
    vint c0(sz);
    cjk = vvint(sz, c0);
    for (int i = 0; i < sz; i++)
    {
        cjk[i][i] = 0;
        for (int j = 0; j < sz; j++)
            if ( i != j )
                cjk[i][j] = modinv(qs[i], qs[j]);
    }
}

rns_ns::RnsMrs::RnsMrs(const RnsMrs & A, Rns::Op op, const RnsMrs & B)
{
    vint myqs = qsetOp(A.getQs(), op, B.getQs());
    Rns::initV(myqs);
    myinit();
}

vint rns_ns::RnsMrs::pow2div(const vint & v, int pow2) const
{
    vint r = v;
    while (pow2--) r = pow2div1(r);
    return r;
}

vint rns_ns::RnsMrs::pow2div1(const vint & v) const
{
    vint m = mrs(v);
    int parity = 0;
    for (auto x : m) parity = (parity + x) % 2;
    RnsForm a(*this, v);
    if (parity) --a;
    auto av = a.values();
    div2exact(av);
    return av;
}

vint rns_ns::RnsMrs::mrs(const vint & x) const
{
    int sz = size();
    vint v(sz), xb(sz);

    for (int k = 0; k < sz; k++)
    {
        auto q = qs[k];
        v[k] = x[k];
        for (int j = 0; j < k; j++)
        {
            auto vkj = modsub(v[k], v[j], q);
            v[k] = modmul(vkj, cjk[j][k], q);
        }
        xb[k] = v[k];
    }
    return xb;
}

void rns_ns::Rns::div2exact(vint & v) const
{
    int sz = size();
    for (int i = 0; i < sz; i++)
    {
        auto & x = v[i];
        if (x % 2) x = (x + qs[i]) / 2;
        else x /= 2;
    }
}

rns_ns::RnsShrinkRound::RnsShrinkRound(const Rns & aQ, const Rns & aP)
    : Q(aQ), P(aP)
{
    // compute P1q, P2q, P2

    // calc P1q
    // plan is to fint P-1 in Q then +1, then invert
    RnsForm Pm1(P, 0), p1(P, 1);
    Pm1 -= p1;
    RnsForm Pm1inQ = Pm1.rebaseAny(Q);
    PinQ = Pm1inQ + RnsForm(Q, 1);
    P1q = PinQ.invert();

    // P2 in P
    P2p = Pm1.div2exact();

    // P2 in Q
    P2q = P2p.rebaseAny(Q);
}
