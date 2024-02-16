#pragma once

#include <string>

#include "integer.h"
#include "err.h"

using std::string;

namespace rns_ns
{

class RnsForm;

class Rns
{
    protected:
        vint qs, Ms_, us;
        Integer maxPow2digit = -1;
        vint neg(vint v) const;
        vint split(Integer x) const;
        vint splitFmod(double x) const;
        Integer lowval(const vint & v) const;
        bool islowval(const vint & v) const;

        std::pair<Integer, int> blend_(const vint & v) const;
        double blendDbl(const vint & v, bool recenter) const;

        void initL(std::initializer_list<Integer> li);
        void initV(const vint & ve);
        Rns(std::initializer_list<Integer> li) { initL(li); }
        Rns() {}
        virtual string print() const;
        virtual vint pow2div(const vint & v, int pow2) const = 0;

    private:
        void assertRnsInited() const { if (qs.empty()) nevers("rns is not set"); }
        friend class RnsForm;
    public:
        int size() const { return (int)qs.size(); }
        Integer dynrange_() const;
        double dynrangeDbl() const;
        Integer q(int i) const { return qs[i]; }
        virtual vint mrs(const vint & v) const = 0;
        vint getQs() const { return qs; }
        void div2exact(vint & v) const;

    public:
        enum Op { minus, plus };
        static vint qsetOp(const vint & a, Op, const vint & b);

};

struct RnsShrinkRound;
class RnsForm
{
        const Rns * prns = nullptr;
        vint v;

        RnsForm(const Rns * r) : prns(r) {}
    public:
        int size() const { return (int)v.size(); }
        RnsForm() {}

        RnsForm(const Rns & r, Integer x) : prns(&r), v(r.split(x)) {}
        RnsForm(const Rns & r, int zero, double val) : prns(&r), v(r.splitFmod(val)) {}
        RnsForm(const Rns & r, const vint & u) : prns(&r), v(u) {}
        RnsForm(const Rns & r, int zero, vint vs);

        RnsForm(const Rns * r, Integer x) : prns(r), v(r->split(x)) {}
        RnsForm(const Rns * r, const vint & u) : prns(r), v(u) {}

        bool match(const RnsForm & b) const;
        bool match(const Rns * pr) const { return prns == pr; }

        vint values() const { return v;  }
        const Rns * rns() const { return prns; }
        Integer lowval() const { if (!prns) never; return prns->lowval(v); }
        bool islowval() const { if (!prns) never; return prns->islowval(v); }
        RnsForm pow2div(int pow2) const
        {
            if (!prns) never;
            RnsForm r(prns);
            r.v = prns->pow2div(v, pow2);
            return r;
        }
        Integer blend_() const { return prns->blend_(v).first; }
        int rank_() const { return prns->blend_(v).second; }
        double blendDbl(bool recenter) const { return prns->blendDbl(v, recenter); }
        RnsForm div2exact() const { auto r = *this; prns->div2exact(r.v); return r; }
        RnsForm invert() const;

        int operator+=(const vint & b); // returns number of oveflows
        RnsForm & operator+=(const rns_ns::RnsForm & b);
        RnsForm operator+(const RnsForm & b) const { RnsForm t = *this;  return t += b; }
        RnsForm & operator*=(const RnsForm & b);
        RnsForm operator*(const RnsForm & b) const { RnsForm t = *this;  return t *= b; }

        //RnsForm& operator/=(const rns_ns::RnsForm& b);
        RnsForm operator/(const RnsForm & b) const { RnsForm q(prns, 0); divABRQ(*this, b, nullptr, &q); return q; }
        RnsForm operator%(const RnsForm & b) const { RnsForm r(prns, 0); divABRQ(*this, b, &r, nullptr); return r; }
        RnsForm operator/(Integer b) const { return *this / RnsForm(prns, b); }
        RnsForm operator%(Integer b) const { return *this % RnsForm(prns, b); }
        void divABRQ(const RnsForm & a, const RnsForm & b, RnsForm * r, RnsForm * q) const;

        RnsForm operator-() const;
        RnsForm & operator-=(const rns_ns::RnsForm & b) { return *this += -b; }
        RnsForm operator-(const rns_ns::RnsForm & b) const { RnsForm t = *this;  return t -= b; }

        RnsForm & operator--() { return *this -= RnsForm(*prns, 1); }
        RnsForm operator--(int) { auto t = *this; --*this; return t; }
        RnsForm & operator++() { return *this += RnsForm(*prns, 1); }
        RnsForm operator++(int) { auto t = *this; ++*this; return t; }

        bool operator==(const RnsForm & b) const { return prns == b.prns && v == b.v; }
        bool operator!=(const RnsForm & b) const { return !(*this == b); }

        RnsForm rebaseAny(const Rns & newRns) const;
        RnsForm rebaseCut(const Rns & newRns) const; // ABCDE -> ABC
        RnsForm rebaseAdd(const Rns & newRns) const // ABC -> ABCDE
        {
            return rebaseAny(newRns); // can be done a bit more efficiently
        }

        // ABCDE -> ABC with round div: Q=ABC P=DE: QP->Q
        ///RnsForm rebaseShrinkRound(const Rns& Q, const Rns& P) const;
        RnsForm rebaseShrinkRound(const RnsShrinkRound & dat) const;
};

struct RnsShrinkRound
{
    const Rns & Q;
    const Rns & P;
    RnsForm P1q, P2q, P2p; // see details in func
    RnsForm PinQ;
    RnsShrinkRound(const Rns & aQ, const Rns & aP);
};

class RnsMrs : public Rns
// ref: (in russian) K. S. Isupov, Vysokoproizvoditel'nye vychislenija ...
{
        using vvint = std::vector<vint>;
        vvint cjk;

        void myinit();
    public:
        RnsMrs() : Rns() {}
        RnsMrs(std::initializer_list<Integer> li) : Rns(li) { myinit(); }
        RnsMrs(const RnsMrs & A, Rns::Op op, const RnsMrs & B);

        void init(std::initializer_list<Integer> li) { Rns::initL(li); myinit(); }

        virtual string print() const override { return Rns::print(); }
        virtual vint pow2div(const vint & v, int pow2) const override;

    private:
        vint pow2div1(const vint & v) const;
        vint mrs(const vint & v) const override;
};

class RnsSel : public Rns
// ref: Selianinau, An Efficient CRT-Base Power-of-Two Scaling ...
{

};

struct RankAlgYesData
{
    ///int n;
    vint v;
    vint xchi, xksi;
    int xrnk_lo, xrnk_hi;
    int xpar_lo, xpar_hi;
    vint ychi, yksi;
    int yrnk, ypar = 0;
    int xrnk, xpar = 0;
    int xrnQ = 0, xpaQ = 0; // opposite to xrnk and xpar
};

class RnsYes : public Rns
// Oleg's algorithm
{
        Integer beta;
        RnsForm Y_, Yn_, Yn1;
        RnsForm Zx_; // alg v5
        vint Zchi; // alg v5
        int rYn1 = -1; // alg v7 rank of Yn1
        vint Yn1chi; // alg v7 chi values of Yn1
        vint qifb, qimb; // qs inverted in beta

        struct V10
        {
            Integer mod, Yn1mod;
            vint mimods;
            Integer Mmod; // M%mod
        } v10 { 0, 0, {}, 0 };

        void myinit();
    public:
        RnsYes() : Rns(), beta(0) {}
        RnsYes(std::initializer_list<Integer> li) : Rns(li) { myinit(); }
        void init(std::initializer_list<Integer> li) { Rns::initL(li); myinit(); }
        virtual string print() const override { return Rns::print() + myprint(); }
        virtual vint pow2div(const vint & v, int pow2) const override;

        RnsForm getY() const { return Y_; }
        vint mrs(const vint & v) const override { never; }

    private:
        vint pow2div1(const vint & v) const;
        void div2exact(vint & v) const;
        int parity_v1(const vint & v, int level) const;
        string myprint() const;
        int parity(const vint & v) const { return parity_v1(v, 2); }

        vint rank_2chi(const vint & x) const;

        vint rank_cube_v2(const vint & chi) const;
        vint rank_cube_v8(const vint & chi) const;
        std::tuple<int, bool, bool> rank_diag_v2(const vint & ksi) const;
        vint rank_cube_round(const vint & chi) const { return rank_cube_v2(chi); }
        std::tuple<int, int, int> rank_diag_round_skew(const vint & ksi) const;
        std::tuple<int, int> rank_diag_round(const vint & ksi) const;
        vint rank_cube_floor(const vint & chi) const;
        std::tuple<int, int> rank_diag_floor(const vint & ksi) const;

        int rank_parity(const vint & chi, int rnk) const;
        vint rank_ypoint(const vint & chi, int k) const;

    public:
        std::pair<int, int> rank_v1(const vint & v) const;
        std::pair<int, int> rank_v2(const vint & v) const;
        std::pair<int, int> rank_v3(const vint & v) const;
        std::pair<int, int> rank_v4(const vint & v) const;
        std::pair<int, int> rank_v5(const vint & v) const;
        std::pair<int, int> rank_v5(RankAlgYesData & data) const;
        void rank_v5_chi(RankAlgYesData & data) const;
        void rank_v5_nodiag(RankAlgYesData & data) const;
        void rank_v5_xrnk(RankAlgYesData & data) const;
        void rank_v5_addz(RankAlgYesData & data) const;
        static int deltaCount(const vint & v);
        vint negate(const vint & v) const { return neg(v); }

        std::pair<int, int> rank_v6(const vint & v) const;
        std::pair<int, int> rank_v6(RankAlgYesData & data) const;

        std::pair<int, int> rank_v7(const vint & v) const;
        int deltaPlus(const vint & a, const vint & b) const;
        std::pair<int, int> rank_v8(const vint & v) const;

        std::pair<int, int> rank_v9(const vint & v) const;
        std::pair<int, int> rank_v9(RankAlgYesData & data) const;
        void rank_v9_chi(RankAlgYesData & data) const { return rank_v5_chi(data); }
        void rank_v9_nodiag(RankAlgYesData & data) const { return rank_v5_nodiag(data); }
        void rank_v9_xrnk(RankAlgYesData & data) const;
        void rank_v9_x2(RankAlgYesData & data) const;

        std::pair<int, int> rank_v10(const vint & v) const;
        Integer rank_rem(const vint & vchi, int rnk) const;

        std::tuple<int, int, int> rank_basic_alg(const vint & v) const;

};

} // rns_ns


///rns_ns::RnsForm operator+(const rns_ns::RnsForm& a, const rns_ns::RnsForm& b);
///rns_ns::RnsForm operator*(const rns_ns::RnsForm& a, const rns_ns::RnsForm& b);
///rns_ns::RnsForm operator-(const rns_ns::RnsForm& a, const rns_ns::RnsForm& b);
///rns_ns::RnsForm operator-(const rns_ns::RnsForm& a);

///inline rns_ns::RnsForm& operator+=(rns_ns::RnsForm& a, const rns_ns::RnsForm& b) { return a = a + b; }
///inline rns_ns::RnsForm& operator*=(rns_ns::RnsForm& a, const rns_ns::RnsForm& b) { return a = a * b; }
///inline rns_ns::RnsForm& operator-=(rns_ns::RnsForm& a, const rns_ns::RnsForm& b) { return a = a - b; }


