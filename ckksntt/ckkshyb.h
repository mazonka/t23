#pragma once

#include "ckkselem.h"

namespace poly
{
using Dpoly = std::vector<poly::Poly>;
} //poly

namespace ckks
{

struct EkHyb
{
    int level;
    Integer P, ql;
    poly::Dpoly db, da;
    EkHyb(int level, Sk sk, Param p, RndStream & rs);
    string print() const;
};

} // ckks


// Polymonial digit decomposition
namespace poly
{
Dpoly WD(const poly::Poly & a, Integer w, Integer q);
Dpoly PW(const poly::Poly & a, Integer w, Integer q);
int calc_dnum(Integer w, Integer q);
poly::Poly dot(const Dpoly & a, const Dpoly & b, Integer q);

} //poly

namespace ckks
{
Ctxt relinHyb(const Ctxt3 & c, const Param & p, const EkHyb & ek);
Ctxt mulHyb(const Ctxt & a, const Ctxt & b, const Param & p, const EkHyb & ek);
} // ckks
