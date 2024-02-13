#include <iostream> /// debug
using std::cout;

#include <sstream>
#include <algorithm>

#include "rns.h"
#include "err.h"


// fixed parity/rank relation
std::pair<int, int> rns_ns::RnsYes::rank_v3(const vint & v) const
{
    int n = size();
    // check if 0
    int delta = deltaCount(v);
    if (0 == delta) return { 0, 0 };


    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_round(xchi);
    auto [xrnk_lo, xrnk_hi] = rank_diag_round(xksi);
    int xpar_lo = rank_parity(xchi, xrnk_lo);
    if (xrnk_lo == xrnk_hi) return { xrnk_lo, xpar_lo };

    vint ychi = rank_ypoint(xchi, n + 1);
    vint yksi = rank_cube_round(ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(yksi);
    if (yrnk_lo != yrnk_hi) never;
    auto yrnk = yrnk_lo;
    int yparity = rank_parity(ychi, yrnk);

    int xpar_hi = rank_parity(xchi, xrnk_hi);
    if (yparity == xpar_lo) return { xrnk_lo, xpar_lo };
    if (yparity == xpar_hi) return { xrnk_hi, xpar_hi };
    never;
}

vint rns_ns::RnsYes::rank_cube_floor(const vint & v) const
{
    const int n = size();
    vint ksi(n);
    for (int i = 0; i < n; i++)
    {
        const bool ORI_MOD_FAST = true;
        if (ORI_MOD_FAST)
        {
            auto m = qs[i];
            auto m1mb = qimb[i];
            // (((a * b) % m) * m1mb) % b;
            auto x1 = modmul(v[i], beta, m);
            auto x5 = modmul(m1mb, x1, beta);
            ksi[i] = x5;

        }
        else
        {
            auto m = qs[i];
            double x1 = v[i] * 1.0;
            double x2 = x1 * beta / m;
            ksi[i] = int(x2 + 0.5);
        }
    }
    return ksi;
}

std::tuple<int, int> rns_ns::RnsYes::rank_diag_floor(const vint & v) const
{
    const int n = size();
    Integer sum = 0;
    int rank = 0;
    for (auto x : v)
    {
        sum += x;
        if (sum >= beta)
        {
            ++rank;
            sum -= beta;
        }
    }
    auto apx = modadd(sum, n, beta);
    bool diag = (apx < n);
    if (!diag) return { rank, rank };
    return { rank, rank + 1 };
}

// replaceed with floor
std::pair<int, int> rns_ns::RnsYes::rank_v4(const vint & v) const
{
    int n = size();

    vint xchi = rank_2chi(v);

    // core 1
    vint xksi = rank_cube_floor(xchi);
    auto [xrnk_lo, xrnk_hi] = rank_diag_floor(xksi);
    int xpar_lo = rank_parity(xchi, xrnk_lo);

    if (xrnk_lo == xrnk_hi) return { xrnk_lo, xpar_lo };
    vint ychi = rank_ypoint(xchi, n + 1);

    // core 2
    vint yksi = rank_cube_floor(ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_floor(yksi);
    int yparity = rank_parity(ychi, yrnk_lo);

    if (yrnk_lo != yrnk_hi) never;

    int xpar_hi = rank_parity(xchi, xrnk_hi);
    if (yparity == xpar_lo) return { xrnk_lo, xpar_lo };
    if (yparity == xpar_hi) return { xrnk_hi, xpar_hi };
    never;
}


std::pair<int, int> rns_ns::RnsYes::rank_v5(const vint & v) const
{
    RankAlgYesData data;
    data.v = v;
    return rank_v5(data);
}

std::pair<int, int> rns_ns::RnsYes::rank_v5(RankAlgYesData & d) const
{
    // compute non-diagonals and exit
    // if diagonals compute Y-offset
    // then check theorem pairing-rho (problematic region)
    // if ok return
    // if not, compute for Z-point
    // Z-point is M2=M/2 offset (outside of problematic region)
    // Note: [-Y,Y] - problematic region
    // Note: -Y+M2 must be > Y; and Y+M2 < M-Y
    // parity of Z detemines if overflow happened

    // check if 0
    int delta = deltaCount(d.v);
    if ( 0 == delta ) return { 0, 0 };

    // compute non-diagonals and exit
    rank_v5_chi(d);
    rank_v5_nodiag(d);
    if (d.xrnk_lo == d.xrnk_hi) return { d.xrnk_lo, d.xpar_lo };

    // ranks are questionable => compute Y point
    rank_v5_xrnk(d);
    //return { d.xrnk, d.xpar }; // this must be the same as v3

    auto u = d;
    u.v = negate(d.v);
    rank_v5_chi(u);
    rank_v5_nodiag(u);
    rank_v5_xrnk(u);

    // check Theorem 1
    if (d.xrnk + u.xrnk == delta - 1)
    {
        // parity must correspond
        if (d.xpar == u.xpar) never;
        return { d.xrnk, d.xpar };
    }

    auto z = d;
    rank_v5_chi(z);
    rank_v5_addz(z);
    rank_v5_nodiag(z);
    if (z.xrnk_lo == z.xrnk_hi)
    {
        // we are off the diagonal
        if (z.xpar_lo == d.xpar_lo) return { d.xrnk_lo, d.xpar_lo };
        if (z.xpar_lo == d.xpar_hi) return { d.xrnk_hi, d.xpar_hi };
        never;
    }

    return { d.xrnk, d.xpar };
    rank_v5_xrnk(z);

    never;
}

void rns_ns::RnsYes::rank_v5_chi(RankAlgYesData & d) const
{
    d.xchi = rank_2chi(d.v);
}

void rns_ns::RnsYes::rank_v5_nodiag(RankAlgYesData & d) const
{
    d.xksi = rank_cube_round(d.xchi);
    std::tie(d.xrnk_lo, d.xrnk_hi) = rank_diag_round(d.xksi);
    // assertions
    {
        if (d.xrnk_hi < d.xrnk_lo) never;
        if (d.xrnk_lo < 0) never;
        if (d.xrnk_hi >= size()) never;
    }
    d.xpar_lo = rank_parity(d.xchi, d.xrnk_lo);
}

void rns_ns::RnsYes::rank_v5_xrnk(RankAlgYesData & d) const
{
    int n = size();
    d.ychi = rank_ypoint(d.xchi, n + 1);
    // why +1
    // parity of Y vector is same as of n because Y=EMi (each Mi is odd)
    // hence parity of Yn1 must be even
    d.yksi = rank_cube_round(d.ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(d.yksi);
    if (yrnk_lo != yrnk_hi) never;
    d.yrnk = yrnk_lo;
    d.ypar = rank_parity(d.ychi, d.yrnk);

    d.xpar_hi = rank_parity(d.xchi, d.xrnk_hi);
    if (0) {}
    else if (d.ypar == d.xpar_lo) { d.xrnk = d.xrnk_lo; d.xpar = d.xpar_lo; }
    else if (d.ypar == d.xpar_hi) { d.xrnk = d.xrnk_hi; d.xpar = d.xpar_hi; }
    else never;
}

void rns_ns::RnsYes::rank_v5_addz(RankAlgYesData & d) const
{
    int n = size();
    for (int i = 0; i < n; i++) d.xchi[i] = modadd(d.xchi[i], Zchi[i], qs[i]);
}

int rns_ns::RnsYes::deltaCount(const vint & v)
{
    int c = 0;
    for ( auto x : v )
        if ( x != 0 ) ++c;

    return c;
}

/*///
vint rns_ns::RnsYes::invert(const vint & v) const
{
    int n = size();
    if (v.size() != n) never;
    vint r(n);
    for (int i = 0; i < n; i++)
        if (v[i] == 0) r[i] = 0;
        else r[i] = qs[i] - v[i];
    return r;
}
*/

std::pair<int, int> rns_ns::RnsYes::rank_v6(const vint & v) const
{
    return std::pair<int, int>();
}

std::pair<int, int> rns_ns::RnsYes::rank_v6(RankAlgYesData & data) const
{
    // multiply by 3 or greater coeff
    return std::pair<int, int>();
}

std::pair<int, int> rns_ns::RnsYes::rank_v7(const vint & v) const
{
    // use Yn1 and theorem 2

    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_round(xchi);
    auto [xrnk_lo, xrnk_hi] = rank_diag_round(xksi);
    int xpar_lo = rank_parity(xchi, xrnk_lo);
    if (xrnk_lo == xrnk_hi) return { xrnk_lo, xpar_lo };

    //cout << " [D] ";

    vint ychi = rank_ypoint(xchi, n + 1);
    vint yksi = rank_cube_round(ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(yksi);
    if (yrnk_lo != yrnk_hi) never;
    auto yrnk = yrnk_lo;

    // Theorem 2 says r(x+y) = r(x)+r(y)-delta(a,b)
    // r(x) = r(x+y) + delta(a,b) - r(y)

    auto rxy = yrnk;
    auto ry = rYn1;
    auto delta = deltaPlus(xchi, Yn1chi);

    auto rx = rxy + delta - ry;
    if (rx == xrnk_lo) return { xrnk_lo, xpar_lo };

    int xpar_hi = rank_parity(xchi, xrnk_hi);
    if (rx == xrnk_hi) return { xrnk_hi, xpar_hi };

    never;
}

int rns_ns::RnsYes::deltaPlus(const vint & a, const vint & b) const
{
    int n = size();
    if (n != (int)a.size()) never;
    if (n != (int)b.size()) never;

    int cntr = 0;
    for (int i = 0; i < n; i++)
        if (a[i] + b[i] >= qs[i]) ++cntr;
    return cntr;
}

std::pair<int, int> rns_ns::RnsYes::rank_v8(const vint & v) const
{
    // use Yn1 and theorem 2 same as v7 but
    // X+Yn1 do not wrap - extended cube

    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_round(xchi);
    auto [xrnk_lo, xrnk_hi] = rank_diag_round(xksi);
    int xpar_lo = rank_parity(xchi, xrnk_lo);
    if (xrnk_lo == xrnk_hi) return { xrnk_lo, xpar_lo };

    vint ychi = rank_ypoint(xchi, n + 1); // no effect, just size
    for (int i = 0; i < n; i++) ychi[i] = xchi[i] + n + 1;

    vint yksi = rank_cube_v8(ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(yksi);
    if (yrnk_lo != yrnk_hi) never;
    auto yrnk = yrnk_lo;

    // Theorem 2 says r(x+y) = r(x)+r(y)-delta(a,b)
    // r(x) = r(x+y) + delta(a,b) - r(y)

    auto rxy = yrnk;
    auto ry = rYn1;
    auto delta = 0 * deltaPlus(xchi, Yn1chi);

    auto rx = rxy + delta - ry;
    int ypar = rank_parity(ychi, yrnk);

    if (rx == xrnk_lo) return { xrnk_lo, xpar_lo };

    int xpar_hi = rank_parity(xchi, xrnk_hi);
    if (rx == xrnk_hi) return { xrnk_hi, xpar_hi };

    never;
}

std::pair<int, int> rns_ns::RnsYes::rank_v9(const vint & v) const
{
    RankAlgYesData data;
    data.v = v;
    return rank_v9(data);
}

std::pair<int, int> rns_ns::RnsYes::rank_v9(RankAlgYesData & d) const
{
    // same as v5
    // but use both Theorems
    // detect overflow by parity of X*2 value

    // check if 0
    int delta = deltaCount(d.v);
    if (0 == delta) return { 0, 0 };

    // compute non-diagonals and exit
    rank_v9_chi(d);
    rank_v9_nodiag(d);
    if (d.xrnk_lo == d.xrnk_hi) return { d.xrnk_lo, d.xpar_lo };

    // ranks are questionable => compute Y point
    rank_v9_xrnk(d);
    //return { d.xrnk, d.xpar }; // this must be the same as v3

    auto u = d;
    u.v = negate(d.v);
    rank_v9_chi(u);
    rank_v9_nodiag(u);
    rank_v9_xrnk(u);

    // check Theorem 1
    if (d.xrnk + u.xrnk == delta - 1)
    {
        never; // just for fun to find this case (should be in 4:13;15;17;19 around 7k) later remove
        // parity must correspond
        if (d.xpar == u.xpar) never;
        return { d.xrnk, d.xpar };
    }

    ///return { d.xrnk, d.xpar };

    auto z = d;
    rank_v9_chi(z);
    rank_v9_x2(z);
    rank_v9_nodiag(z);
    if (z.xrnk_lo != z.xrnk_hi)
    {
        // we are on the diagonal
        rank_v9_xrnk(z);
        ///if (z.xpar == 0) return { d.xrnk, d.xpar };
        ///else return { d.xrnQ, d.xpaQ };
    }
    else
    {
        z.xrnk = z.xrnk_lo; z.xpar = z.xpar_lo;
    }

    //never;
    //return { d.xrnk, d.xpar };

    // now check T2
    // if x+x<M: r(z) = r(x+x) = r(x)+r(x)-delta(x,x)

    auto rz = z.xrnk;
    auto rx = d.xrnk;
    auto delta2 = deltaPlus(d.xchi, d.xchi);

    auto r2 = rx * 2 - delta2;

    if (rz == r2) // no overflow
        return { d.xrnk, d.xpar };
    else // overflow
        return { d.xrnQ, d.xpaQ };
}

void rns_ns::RnsYes::rank_v9_xrnk(RankAlgYesData & d) const
{
    int n = size();
    d.ychi = rank_ypoint(d.xchi, n + 1);
    // why +1
    // parity of Y vector is same as of n because Y=EMi (each Mi is odd)
    // hence parity of Yn1 must be even
    d.yksi = rank_cube_round(d.ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(d.yksi);
    if (yrnk_lo != yrnk_hi) never;
    d.yrnk = yrnk_lo;

    // Theorem 2: if x+y<M: r(x+y) = r(x)+r(y)-delta(x,y)
    // r(x) = r(x+y) + delta(x,y) - r(y)

    auto rxy = d.yrnk;
    auto ry = rYn1;
    auto delta = deltaPlus(d.xchi, Yn1chi);

    auto rx = rxy + delta - ry;
    d.xpar_hi = rank_parity(d.xchi, d.xrnk_hi);
    d.ypar = rank_parity(d.ychi, d.yrnk); // we do not need if using T2

    if (0) {}
    else if (rx == d.xrnk_lo)
    {
        d.xrnk = d.xrnk_lo; d.xpar = d.xpar_lo;
        d.xrnQ = d.xrnk_hi; d.xpaQ = d.xpar_hi;
    }
    else if (rx == d.xrnk_hi)
    {
        d.xrnk = d.xrnk_hi; d.xpar = d.xpar_hi;
        d.xrnQ = d.xrnk_lo; d.xpaQ = d.xpar_lo;
    }
    else never;
}

void rns_ns::RnsYes::rank_v9_x2(RankAlgYesData & d) const
{
    int n = size();
    for (int i = 0; i < n; i++) d.xchi[i] = modadd(d.xchi[i], d.xchi[i], qs[i]);
}

std::pair<int, int> rns_ns::RnsYes::rank_v10(const vint & v) const
{
    // use Yn1 and remainder of a coprime test number

    // check if 0
    ///int delta = deltaCount(v);
    if (0 == deltaCount(v)) return { 0, 0 };

    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_round(xchi);
    auto [xrnk_lo, xrnk_hi] = rank_diag_round(xksi);
    int xpar_lo = rank_parity(xchi, xrnk_lo);
    if (xrnk_lo == xrnk_hi) return { xrnk_lo, xpar_lo };

    //cout << " [D] ";

    vint ychi = rank_ypoint(xchi, n + 1);
    vint yksi = rank_cube_round(ychi);
    auto [yrnk_lo, yrnk_hi] = rank_diag_round(yksi);
    if (yrnk_lo != yrnk_hi) never;
    auto yrnk = yrnk_lo;

    // compute remainder for xrnk_lo and xrnk_hi
    auto xrem_lo = rank_rem(xchi, xrnk_lo);
    auto xrem_hi = rank_rem(xchi, xrnk_hi); // can be optimized on xchi
    auto yrem = rank_rem(ychi, yrnk);

    auto xadd_lo = modadd(xrem_lo, v10.Yn1mod, v10.mod);
    auto xadd_hi = modadd(xrem_hi, v10.Yn1mod, v10.mod);

    if ( xadd_lo == yrem ) return { xrnk_lo, xpar_lo };

    int xpar_hi = rank_parity(xchi, xrnk_hi);
    if ( xadd_hi == yrem ) return { xrnk_hi, xpar_hi };

    never;
}

Integer rns_ns::RnsYes::rank_rem(const vint & vchi, int rnk) const
{
    Integer sum { 0 };

    int n = size();
    if (n != (int)v10.mimods.size()) never;

    for (int i = 0; i < n; i++)
    {
        auto t = modmul(v10.mimods[i], vchi[i], v10.mod);
        sum = modadd(t, sum, v10.mod);
    }

    auto m = modmul(v10.Mmod, Integer { rnk }, v10.mod);
    sum = modsub(sum, m, v10.mod);
    return sum;
}

std::tuple<int, int, int> rns_ns::RnsYes::rank_basic_alg(const vint & v) const
{
    if (0 == deltaCount(v)) return { 0, 0, 0 };

    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_round(xchi);
    auto ret = rank_diag_round_skew(xksi);
    return ret;
}

void rns_ns::RnsYes::myinit()
{
    //beta;
    beta = 1;
    Integer m = qs.back();
    while (!!m) { m >>= 1; beta <<= 1; }

    //RnsForm Y, Y2;
    int sz = size();
    RnsForm y(this, 0);

    for (auto q : qs)
    {
        RnsForm m(*this, 1);
        for (auto s : qs)
            if (s != q)
            {
                m *= RnsForm(this, s);
            }
        y += m;
    }
    Y_ = y;
    Yn_ = y * RnsForm(this, sz);
    Yn1 = Yn_ + Y_;

    //vint qib; // qs inverted in beta
    for (int i = 0; i < sz; i++)
    {
        auto mi = modinv(qs[i], beta);
        qifb.push_back(mi);
        qimb.push_back(beta - mi);
    }

    // alg v5
    // compute Z ~ M/2 and even
    Integer q4 = 1;
    for (auto q : qs)
    {
        q4 *= q % Integer(4);
        q4 %= Integer(4);
    }
    RnsForm Za(this, qs);
    Za -= RnsForm(this, q4);
    vint Zv = Za.values();
    div2exact(Zv);
    Zx_ = RnsForm(this, Zv);
    Zchi = rank_2chi(Zv);

    // alg v7
    auto Yn1v = Yn1.values();
    rYn1 = rank_v7(Yn1v).first;
    Yn1chi = rank_2chi(Yn1v);

    // alg v9
    // validate condition Yn1<M/4
    const bool warn = true;
    auto minq = 4 * sz * (sz + 1);
    if (qs[0] < minq)
    {
        string err = "condition 9 (Yn1<M/4) is not satistied";
        if (warn)
        {
            cout << "Warning: " << err << '\n';
            cout << "         n=" << sz << " and q[0]="
                 << qs[0] << " but 4n(n+1)=" << minq << '\n';
        }
        else nevers(err);
    }

    // alg 10 init
    const bool forceV10 = !true;
    v10.mod = Integer { 3 };
    while (1)
    {
        bool ok = true;
        for (auto q : qs)
        {
            if ( !coprime(q, v10.mod) )
            {
                ok = false;
                break;
            }
        }
        if (ok) break;
        v10.mod += 2;
        if (v10.mod > qs.back())
        {
            string err = "cannot find v10 test modulus";
            if ( forceV10 ) nevers(err);
            cout << "Warning: " << err << '\n';
            return;
        }
    }

    v10.Mmod = 1;
    for (auto q : qs)
    {
        auto qm = q % v10.mod;
        v10.Mmod = modmul(v10.Mmod, qm, v10.mod);
        Integer mi = 1;
        for (auto p : qs)
        {
            if (p == q) continue;
            mi = modmul(mi, p, v10.mod);
        }
        v10.mimods.push_back(mi);
    }

    // now its safe to get rem Yn1
    v10.Yn1mod = rank_rem(Yn1chi, rYn1);
}

vint rns_ns::RnsYes::pow2div(const vint & v, int pow2) const
{
    vint r = v;
    while (pow2--) r = pow2div1(r);
    return r;
}

vint rns_ns::RnsYes::pow2div1(const vint & v) const
{
    int par = parity(v);
    RnsForm a(*this, v);
    if (par) --a;
    auto av = a.values();
    div2exact(av);
    return av;
}

void rns_ns::RnsYes::div2exact(vint & v) const
{
    int sz = size();
    for (int i = 0; i < sz; i++)
    {
        auto & x = v[i];
        if (x % 2) x = (x + qs[i]) / 2;
        else x /= 2;
    }
}

int rns_ns::RnsYes::parity_v1(const vint & v, int level) const
{
    if (level == 0) never;
    int sz = size();

    // step 1: compute normalized coordinates
    vint xi(sz);
    for (int i = 0; i < sz; i++) xi[i] = modmul(v[i], us[i], qs[i]);

    // step 2: compute extended coordinates
    vint e(sz);
    for (int i = 0; i < sz; i++)
    {
        auto m = qs[i];
        auto m2 = m >> 1;
        auto m1b = qifb[i];
        // (m1b * (m2 + b - ((a * b + m2) % m)) % b) % b
        auto x1 = modmul(xi[i], beta, m);
        auto x2 = modadd(x1, m2, m);
        auto x3 = beta - x2;
        auto x4 = modadd(m2, x3, beta);
        auto x5 = modmul(m1b, x4, beta);
        e[i] = x5;
    }

    // step 3: check if value is diagonal
    Integer esumMod = 0, esumInt = 0;
    int nb = 0;
    for (auto x : e)
    {
        esumMod = modadd(esumMod, x, beta);
        esumInt += x;
        if (esumInt >= beta) { ++nb; esumInt -= beta; }
    }
    auto apx = modadd(esumMod, sz / 2, beta);
    bool diagonal = (apx < sz);

    if (diagonal)
    {
        RnsForm Xp(*this, v); // use xi instead of x;
        // addition is the same in norm-coordinates
        Xp += Y_; // we might need add Yn=Y*sz and the range Dyn-Yn ?
        auto r = parity_v1(Xp.values(), level - 1);
        if (sz % 2)
            return 1 - r;
        return r;
    }

    // step 4: compute parity for non-diagonal element
    int par = 0;
    for (auto x : xi) par = (par + x) % 2;
    par += nb;

    return par % 2;
}

string rns_ns::RnsYes::myprint() const
{
    std::ostringstream os;
    os << "beta = " << beta << '\n';
    os << "Y ="; for (auto x : Y_.values()) os << ' ' << x; os << '\n';
    os << "Yn ="; for (auto x : Yn_.values()) os << ' ' << x; os << '\n';
    os << "Yn1 ="; for (auto x : Yn1.values()) os << ' ' << x; os << '\n';
    os << "Z ="; for (auto x : Zx_.values()) os << ' ' << x; os << '\n';
    os << "Zchi ="; for (auto x : Zchi) os << ' ' << x; os << '\n';
    os << "Yn1chi ="; for (auto x : Yn1chi) os << ' ' << x; os << '\n';
    os << "qifb ="; for (auto x : qifb) os << ' ' << x; os << '\n';
    os << "qimb ="; for (auto x : qimb) os << ' ' << x; os << '\n';

    os << "Y,Yn = " << blend_(Y_.values()).first;
    os << ' ' << blend_(Yn_.values()).first; os << '\n';

    os << "Yn1,Z = " << blend_(Yn1.values()).first
       << " (rank=" << rYn1 << "), ";
    os << ' ' << blend_(Zx_.values()).first; os << '\n';

    return os.str();
}

// Yes algorithm
std::pair<int, int> rns_ns::RnsYes::rank_v1(const vint & v) const
{
    int n = size();
    Integer iYn { n };

    // step 1: compute normalized coordinates
    vint xchi(n);
    for (int i = 0; i < n; i++) xchi[i] = modmul(v[i], us[i], qs[i]);

    // step 2: compute extended coordinates
    auto extend2cube = [n, this](const vint & chi) -> vint
    {
        vint ksi(n);
        for (int i = 0; i < n; i++)
        {
            const bool ORI_MOD_FAST = true;
            if (ORI_MOD_FAST)
            {
                auto m = qs[i];
                auto m2 = m >> 1;
                auto m1b = qifb[i];
                // (m1b * (m2 + b - ((a * b + m2) % m)) % b) % b
                auto x1 = modmul(chi[i], beta, m);
                auto x2 = modadd(x1, m2, m);
                auto x3 = beta - x2;
                auto x4 = modadd(m2, x3, beta);
                auto x5 = modmul(m1b, x4, beta);
                ksi[i] = x5;

            }
            else
            {
                auto m = qs[i];
                double x1 = chi[i] * 1.0;
                double x2 = x1 * beta / m;
                ksi[i] = int(x2 + 0.5);
            }
        }
        return ksi;
    };
    vint xksi = extend2cube(xchi);

    // step 3: estimate rank and check if value is diagonal
    Integer esumModx = 0, esumIntx = 0;
    int assumed_rankx = 0;
    for (auto x : xksi)
    {
        esumModx = modadd(esumModx, x, beta);
        esumIntx += x;
        if (esumIntx >= beta)
        {
            ++assumed_rankx;
            esumIntx -= beta;
        }
    }
    auto apxx = modadd(esumModx, n / 2, beta);
    bool diagonalx = (apxx < n);

    // step 4: compute parity for non-diagonal element
    int assumed_xpar = 0;
    for (auto x : xchi) assumed_xpar = (assumed_xpar + x) % 2;
    assumed_xpar += assumed_rankx;
    assumed_xpar %= 2;

    if (!diagonalx) return { assumed_rankx, assumed_xpar };

    bool above = (esumIntx < beta / 2);

    // step 5: compute Y-point and its overflows
    ///*RnsForm ypnt = Yn;
    ///*int ovfs = (ypnt += v);

    // step 6: compute chi of Y-point
    vint ychi(n);
    for (int i = 0; i < n; i++) ychi[i] = modadd(xchi[i], iYn, qs[i]);
    //for (int i = 0; i < n; i++) ychi[i] = xchi[i] + iYn;

    if (0) // DBG check that ychi correspond to ypnt
    {
        RnsForm ypnt = Yn_;
        int ovfs = (ypnt += v);
        vint tchi(n);
        vint ypntv = ypnt.values();
        for (int i = 0; i < n; i++) tchi[i] = modmul(ypntv[i], us[i], qs[i]);
        for (int i = 0; i < n; i++) if (tchi[i] != ychi[i]) never;
    }

    // step 7: compute ksi for Y-point
    vint yksi = extend2cube(ychi);

    // step 8: compute parity of Y-point
    Integer esumMody = 0, esumInty = 0; // later: replace esumInty with debug check
    int y_rank = 0;
    for (auto x : yksi)
    {
        esumMody = modadd(esumMody, x, beta);
        esumInty += x;
        if (esumInty >= beta)
        {
            ++y_rank;
            esumInty -= beta;
        }
    }
    auto apxy = modadd(esumMody, n / 2, beta);
    bool diagonaly = (apxy < n);

    if (diagonaly) never;

    int ypar = 0;
    for (auto x : ychi) ypar = (ypar + x) % 2;
    ypar += y_rank;
    //ypar += ovfs;
    ypar %= 2;

    // step 9: determine rankof X-point
    // parity of Y vector is same as of n because Y=EM (each M is odd)
    // parity of Yn is same as Y
    int xpar = (ypar + n) % 2;

    if (xpar == assumed_xpar)
    {
        if (y_rank == assumed_rankx)
            return { assumed_rankx, xpar };
        else
            return { assumed_rankx, xpar };
    }
    else
    {
        if (y_rank == assumed_rankx)
            return { assumed_rankx, assumed_xpar };
        else
        {} // never;
    }

    int rank = assumed_rankx;
    if (above) --rank;
    else ++rank;

    return { rank, xpar };
}


vint rns_ns::RnsYes::rank_2chi(const vint & v) const
{
    const int n = size();
    vint chi(n);
    for (int i = 0; i < n; i++) chi[i] = modmul(v[i], us[i], qs[i]);
    return chi;
}

vint rns_ns::RnsYes::rank_cube_v2(const vint & v) const
{
    const int n = size();
    vint ksi(n);
    for (int i = 0; i < n; i++)
    {
        const bool ORI_MOD_FAST = true;
        if (ORI_MOD_FAST)
        {
            auto m = qs[i];
            auto m2 = m >> 1;
            auto m1b = qifb[i];
            // (m1b * (m2 + b - ((a * b + m2) % m)) % b) % b
            auto x1 = modmul(v[i], beta, m);
            auto x2 = modadd(x1, m2, m);
            auto x3 = beta - x2;
            auto x4 = modadd(m2, x3, beta);
            auto x5 = modmul(m1b, x4, beta);
            ksi[i] = x5;

        }
        else
        {
            auto m = qs[i];
            double x1 = v[i] * 1.0;
            double x2 = x1 * beta / m;
            ksi[i] = int(x2 + 0.5);
        }
    }
    return ksi;
}

// extended cube
vint rns_ns::RnsYes::rank_cube_v8(const vint & v) const
{
    const int n = size();
    vint ksi(n);
    for (int i = 0; i < n; i++)
    {
        const bool ORI_MOD_FAST = !true;
        if (ORI_MOD_FAST)
        {
            never; // not implemented
            auto m = qs[i];
            auto m2 = m >> 1;
            auto m1b = qifb[i];
            // (m1b * (m2 + b - ((a * b + m2) % m)) % b) % b
            auto x1 = modmul(v[i], beta, m);
            auto x2 = modadd(x1, m2, m);
            auto x3 = beta - x2;
            auto x4 = modadd(m2, x3, beta);
            auto x5 = modmul(m1b, x4, beta);
            ksi[i] = x5;

        }
        else
        {
            auto m = qs[i];
            double x1 = v[i] * 1.0;
            double x2 = x1 * beta / m;
            ksi[i] = int(x2 + 0.5);
        }
    }
    return ksi;
}

std::tuple<int, bool, bool> rns_ns::RnsYes::rank_diag_v2(const vint & v) const
{
    const int n = size();
    Integer sum = 0;
    int rank = 0;
    for (auto x : v)
    {
        sum += x;
        if (sum >= beta)
        {
            ++rank;
            sum -= beta;
        }
    }
    auto apx = modadd(sum, n / 2, beta);
    bool diag = (apx < n);
    bool above = (sum < (beta / 2));
    return { rank, diag, above };
}

std::tuple<int, int> rns_ns::RnsYes::rank_diag_round(const vint & v) const
{
    auto [a, b, c] = rank_diag_round_skew(v);
    return { a, b };
}

std::tuple<int, int, int> rns_ns::RnsYes::rank_diag_round_skew(const vint & v) const
{
    const int n = size();
    Integer sum = 0;
    int rank = 0;
    for (auto x : v)
    {
        sum += x;
        if (sum >= beta)
        {
            ++rank;
            sum -= beta;
        }
    }
    int n2 = n / 2;
    auto apx = modadd(sum, n2, beta);
    bool diag = (apx < n);
    if (!diag) return { rank, rank, 0 };
    bool aboveDiagonal = (sum < (beta / 2));
    int iapx = int(apx);
    iapx -= n2;
    if (aboveDiagonal) return { rank - 1, rank, iapx }; // can be negative for 0
    return { rank, rank + 1, iapx };
}

int rns_ns::RnsYes::rank_parity(const vint & chi, int rnk) const
{
    int par = 0;
    for (auto x : chi) par = (par + x) % 2;
    par += rnk;
    par %= 2;
    return par;
}

vint rns_ns::RnsYes::rank_ypoint(const vint & chi, int k) const
{
    int n = size();
    Integer iYn1 {k};
    vint ychi(n);
    for (int i = 0; i < n; i++) ychi[i] = modadd(chi[i], iYn1, qs[i]);
    return ychi;
}

std::pair<int, int> rns_ns::RnsYes::rank_v2(const vint & v) const
{
    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube_v2(xchi);
    auto [xrnk, xdiag, above] = rank_diag_v2(xksi);
    int xparity = rank_parity(xchi, xrnk);
    if (!xdiag) return { xrnk, xparity };

    vint ychi = rank_ypoint(xchi, n + 1);
    vint yksi = rank_cube_v2(ychi);
    auto [yrnk, ydiag, yabove] = rank_diag_v2(yksi);
    int yparity = rank_parity(ychi, yrnk);

    if (xparity == yparity)
        return { xrnk, xparity };

    if (above)
        return { xrnk - 1, xparity };

    return { xrnk + 1, xparity };
}
