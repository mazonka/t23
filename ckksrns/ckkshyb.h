#pragma once

#include "ckkselem.h"

namespace poly
{
using Dpoly = std::vector<poly::Poly>;
} //poly

namespace ckks
{

struct EkHybP
{
    int level;
    Integer P, ql;
    poly::Dpoly db, da;
    EkHybP(int level, SkP sk, Param p, RndStream & rs);
    string print() const;
};

struct EkHybR
{
    ///int level;
    Integer P, ql;
    poly::Dpoly db, da;
    EkHybR(SkR sk, Param p, RndStream & rs, rns_ns::Rns & rext, rns_ns::RnsShrinkRound rshrink);
    string print() const;

    static Integer findExtDigit(const vint & v, int n);
};

} // ckks


// Polymonial digit decomposition
namespace poly
{
Dpoly WDp(const poly::Poly & a, Integer w, Integer q);
Dpoly PWp(const poly::Poly & a, Integer w, Integer q);
int calc_dnumP(Integer w, Integer q);
Poly dotP(const Dpoly & a, const Dpoly & b, Integer q);

PolyRns WDr(const PolyRns & a);
PolyRns PWr(const PolyRns & a);
PolyRns dotR(const PolyRns & a, const PolyRns & b);

} //poly

namespace ckks
{
CtxtP relinHybP(const Ctxt3P & c, const Param & p, const EkHybP & ek);
CtxtP mulHybP(const CtxtP & a, const CtxtP & b, const Param & p, const EkHybP & ek);
} // ckks

