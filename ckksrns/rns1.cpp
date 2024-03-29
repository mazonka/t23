#include <iostream> /// debug
using std::cout;

#include <sstream>
#include <algorithm>

#include "rns.h"
#include "err.h"


vint rns_ns::Rns::split(Integer x) const
{
    assertRnsInited();
    vint r;
    for (auto q : qs) r.push_back(x % q);
    return r;
}

Integer rns_ns::Rns::lowval(const vint & v) const
{
    if (!islowval(v)) nevers("rns overflow: cannot represent in q0 range");
    return v[size() - 1];
}

bool rns_ns::Rns::islowval(const vint & v) const
{
    int sz = (int)v.size();
    for (int i = 1; i < sz; i++) // 1
        if (v[0] != v[i]) goto chkmax;
    return true;

chkmax:
    int sz1 = sz - 1;
    for (int i = 0; i < sz1; i++)
        if (v[sz1] % qs[i] != v[i] ) return false;
    return true;
}

void rns_ns::Rns::init(std::initializer_list<Integer> li)
{
    qs = li;
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

rns_ns::RnsForm rns_ns::RnsForm::operator-() const
{
    if (!prns) never;
    const auto & qs = prns->qs;
    int sz = prns->size();
    RnsForm r(*this);
    for (int i = 0; i < sz; i++) r.v[i] = qs[i] - v[i];
    return r;
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
        for (int j = 0; j < sz; j++) // FIXME j<k
        {
            if (k > j)
            {
                auto vkj = modsub(v[k], v[j], q);
                v[k] = modmul(vkj, cjk[j][k], q);
            }
        }
        xb[k] = v[k];
    }
    return xb;
}

void rns_ns::RnsMrs::div2exact(vint & v) const
{
    int sz = size();
    for (int i = 0; i < sz; i++)
    {
        auto & x = v[i];
        if (x % 2) x = (x + qs[i]) / 2;
        else x /= 2;
    }
}

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

    //vint qib; // qs inverted in beta
    for (int i = 0; i < sz; i++)
    {
        auto mi = modinv(qs[i], beta);
        qib.push_back(mi);
    }
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
        auto m1b = qib[i];
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
        RnsForm Xp(*this, v); // FIXME use xi instead of x; addition is the same in norm-coordinates
        Xp += Y_; // we might need add Yn=Y*sz and the range ?? Dyn-Yn ?
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
    os << "qib ="; for (auto x : qib) os << ' ' << x; os << '\n';

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
                auto m1b = qib[i];
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
    Integer esumMody = 0, esumInty = 0; // FIXME remove esumInty - replace with debug check
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

std::pair<int, int> rns_ns::RnsYes::rank_v2(const vint & v) const
{
    int n = size();

    vint xchi = rank_2chi(v);
    vint xksi = rank_cube(xchi);
    auto [xrnk, xdiag, abogit ve] = rank_diag(xksi);
    int xparity = rank_parity(xchi, xrnk);
    if (!xdiag) return { xrnk, xparity };

    vint ychi = rank_ypoint(xchi);
    vint yksi = rank_cube(ychi);
    auto [yrnk, ydiag, yabove] = rank_diag(yksi);
    int yparity = rank_parity(ychi, yrnk);

    if (xparity == yparity)
        return { xrnk, xparity };

    if ( above )
        return { xrnk - 1, xparity };

    return { xrnk + 1, xparity };
}

vint rns_ns::RnsYes::rank_2chi(const vint & v) const
{
    const int n = size();
    vint chi(n);
    for (int i = 0; i < n; i++) chi[i] = modmul(v[i], us[i], qs[i]);
    return chi;
}

vint rns_ns::RnsYes::rank_cube(const vint & v) const
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
            auto m1b = qib[i];
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

std::tuple<int, int> rns_ns::RnsYes::rank_diag(const vint & v) const
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
    if (!diag) return { rank, rank };
    bool above = (sum < (beta / 2));
    if ( above ) return { rank - 1, rank };
    return { rank, rank + 1 };
}

int rns_ns::RnsYes::rank_parity(const vint & chi, int rnk) const
{
    int par = 0;
    for (auto x : chi) par = (par + x) % 2;
    par += rnk;
    par %= 2;
    return par;
}

vint rns_ns::RnsYes::rank_ypoint(const vint & chi) const
{
    int n = size();
    Integer iYn1 { n + 1 };
    vint ychi(n);
    for (int i = 0; i < n; i++) ychi[i] = modadd(chi[i], iYn1, qs[i]);
    return ychi;
}
