#pragma once

#include <map>

#include "bigun.h"
#include "poly.h"


namespace ntt
{
extern bool disabled;


Integer findCloseQ(int n, Integer x);

poly::Poly nttBfly(poly::Poly a, Integer q);
poly::Poly ittBfly(poly::Poly p, Integer q);

///Num find1NthRoot(int n, bool throws = true);
Integer find2NthRoot(Integer mod, int n, bool throws = true);
bool chk1NthRoot(Integer mod, int n, Integer om);
bool chk2NthRoot(Integer mod, int n, Integer ps);

class Context
{
    protected:
        int n = 0;
        Integer mod;
        Integer npsi;
        ///Context1ntt ctx1;
        std::vector<Integer> psi_powers, ps1_powers;
        int d = -1;
        Integer n1;

        void build_psi_powers();
        void calc_d() { d = calcPow2(n); }
        void calc_n1() { n1 = mod::inv(Integer(n), mod); }
        static Integer find_psi(Integer mod, int n);
        void init();

    public:
        Context() : mod(0), npsi(0), n1(0) {}
        Context(int n, Integer m, Integer apsi);
        Context(int n, Integer m);

        Integer psi(int i) const { return psi_powers[i]; }
        Integer ps1(int i) const { return ps1_powers[i]; }
        Integer psi_rev(int i) const { return psi_powers[rev(i)]; }
        Integer ps1_rev(int i) const { return ps1_powers[rev(i)]; }

        //...
        static int calcPow2(int i); // -1 if not
        int rev(int x) const;

        //...
        void assertPow2() const { if (d < 0) throw "n is not power of 2"; }
        //...
        Integer getN1() const { return n1; }
};

class NttMan
{
        using mic = std::map<Integer, Context>;
        using m2m = std::map<int, mic>;
        m2m m;
    public:
        const Context & give(int n, Integer q);
};

} // ntt
