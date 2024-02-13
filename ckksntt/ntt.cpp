
#include "ntt.h"
///#include "ckkselem.h"
#include "err.h"
#include "poly.h"
#include "mathut.h"

namespace g
{
ntt::NttMan nttman;
} // g

Integer ntt::findCloseQ(int n, Integer x)
{
    Integer i(0), U(1), in(n);
    Integer max = x / 2;
    for (; i < max; ++i)
    {
        int two = (i == 0 ? 1 : 2);
        for (int j = 0; j < two; ++j)
        {
            auto mod = (j == 0 ? x + i : x - i);
            if ( math::gcdT(mod, in) != U) continue;
            Integer y = find2NthRoot(mod, n, false);
            if (!!y) return mod;
        }
    }
    nevers("Cannot find close Q for " + std::to_string(todbl(x)));
}

poly::Poly ntt::nttBfly(poly::Poly a, Integer mod)
{
    int n = a.size();
    Context ct = g::nttman.give(n, mod);
    int t = n;

    ct.assertPow2();

    for (int m = 1; m < n; m *= 2)
    {
        t = t / 2;
        for (int i = 0; i < m; i++)
        {
            int j1 = 2 * i * t;
            int j2 = j1 + t - 1;
            auto S = ct.psi_rev(m + i);
            for (int j = j1; j <= j2; ++j)
            {
                auto U = a.v[j];
                //auto V = a.v[j + t] * S;
                //a.v[j] = U + V;
                //a.v[j + t] = U - V;
                auto V = mod::mul(a.v[j + t], S, mod);
                a.v[j] = mod::add(U, V, mod);
                a.v[j + t] = mod::sub(U, V, mod);
            }
        }
    }

    return a;
}

poly::Poly ntt::ittBfly(poly::Poly a, Integer mod)
{
    int n = a.size();
    Context ct = g::nttman.give(n, mod);

    ct.assertPow2();

    int t = 1;
    for (int m = n; m > 1; m /= 2)
    {
        int j1 = 0;
        int h = m / 2;
        for ( int i = 0; i < h; ++i )
        {
            int j2 = j1 + t - 1;
            auto S = ct.ps1_rev(h + i);
            for ( int j = j1; j <= j2; j++ )
            {
                auto U = a.v[j];
                auto V = a.v[j + t];
                //a.v[j] = U + V;
                a.v[j] = mod::add(U, V, mod);
                //a.v[j + t] = (U - V) * S;
                a.v[j + t] = mod::mul( mod::sub(U, V, mod), S, mod);
            }
            j1 = j1 + 2 * t;
        }
        t = 2 * t;
    }

    auto n1 = ct.getN1();

    for ( int j = 0; j < n; j++ )
        //a.v[j] = a.v[j] * n1;
        a.v[j] = mod::mul(a.v[j], n1, mod);

    return a;
}

Integer ntt::find2NthRoot(Integer mod, int n, bool throws)
{
    const auto & m = mod;
    if (m > Integer(1))
    {
        Integer i(2), u(1), mu(mod - u);
        for (; i < m; ++i)
        {
            auto x = mod::pow(i, Integer(n), m); ///i.powmod(n);
            if (x == mu)
            {
                auto om = mod::mul(i, i, m); ///i * i; // .mulmod(i, m);
                if (om == u) continue;
                if (!chk1NthRoot(mod, n, om)) continue;
                return i;
            }
        }
    }

    if (throws)
        nevers("Cannot find root for n=" + std::to_string(n)
               + " mod=" + std::to_string(todbl(mod)));

    return Integer(0);

}

bool ntt::chk1NthRoot(Integer mod, int n, Integer om)
{
    if (om < mod) {}
    else return false;
    Integer u(1);

    if (mod::pow(om, Integer(n), mod) != u) return false;

    if (n % 2 == 0)
    {
        if (mod::pow(om, Integer(n / 2), mod) == (mod - 1)) return true;
    }

    // check that all powers <n !=1
    for (int j = 1; j < n; ++j)
    {
        if (mod::pow(om, Integer(j), mod) == u) return false;
    }

    return true;

}

bool ntt::chk2NthRoot(Integer mod, int n, Integer psi)
{
    Integer i(2), u(1), mu(mod - u);

    auto x = mod::pow(psi, Integer(n), mod);
    if (x == mu)
    {
        auto om = mod::mul(psi, psi, mod);
        if (om == u) return false;
        if (!chk1NthRoot(mod, n, om)) return false;
        return true;
    }

    return false;

}

const ntt::Context & ntt::NttMan::give(int n, Integer q)
{
    auto it1 = m.find(n);
    if (it1 == m.end()) goto build;

    {
        mic & mm = it1->second;
        auto it2 = mm.find(q);
        if (it2 == mm.end()) goto build;
        return it2->second;
    }

build:
    m[n][q] = Context(n, q);
    return m[n][q];
}

void ntt::Context::build_psi_powers()
{
    psi_powers.clear();
    ps1_powers.clear();

    for (int i = 0; i < n; i++)
        psi_powers.push_back(mod::pow(npsi, Integer(i), mod));
    ///npsi.powmod(i));

    auto psi1 = mod::inv(npsi, mod); ///npsi.invmod();
    for (int i = 0; i < n; i++)
        ps1_powers.push_back(mod::pow(psi1, Integer(i), mod));
    ///psi1.powmod(i));
}

Integer ntt::Context::find_psi(Integer mod, int n)
{
    if (n < 1) never;
    if (n == 1) return Integer(1);
    Integer apsi = find2NthRoot(mod, n);
    if (!chk2NthRoot(mod, n, apsi)) nevers("bad psi " + std::to_string(todbl(apsi)) );
    return apsi;
}

void ntt::Context::init()
{
    build_psi_powers();
    calc_d();
    calc_n1();
}

ntt::Context::Context(int an, Integer am, Integer apsi)
    : n(an), mod(am), npsi(apsi), n1(0)///, ctx1()
{
    if (!npsi) return;
    //ctx1 = Context1ntt(an, npsi * npsi);
    init();
}

ntt::Context::Context(int an, Integer am)
    : n(an), mod(am), npsi(find_psi(am, an)), n1(0) ///, ctx1()
{
    //ctx1 = Context1ntt(an, npsi * npsi);
    init();
}

int ntt::Context::calcPow2(int i)
{
    int j = 0;
    for (; i > 1; i /= 2, ++j)
        if (i % 2) return -1;

    return j;
}

int ntt::Context::rev(int x) const
{
    int rv = 0;

    for (int i = 0; i < d; ++i, x >>= 1)
    {
        rv = (rv << 1) | (x & 1);
    }

    return rv;
}
