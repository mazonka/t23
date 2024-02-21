#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>

#include "ckkselem.h"
#include "ntt.h"
#include "chron.h"
#include "main.h"
#include "ckkshyb.h"


using std::string;
using std::cout;
using std::vector;

using poly::Poly;
using poly::PolyRns;



void t05_mul2()
{
    // switch off ntt to
    ntt::Nttoff nttoff; // FIXME find why

    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    using rns_ns::RnsMrs;

    Integer delta_(1024);
    Param param(4, Integer(1024), Integer(delta_), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 0.8, 0.5 };
    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 1.5625, 1.5625 };

    cout << "a =" << a << '\n';

    RnsMrs rns(param.vqs);
    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    auto c3prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << c.c2 << '\n';
    };
    auto c2prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << '\n';
    };

    CtxtP cap = encryptP(skp, map, param, rsP);
    CtxtR car = encryptR(skr, mar, param, rsR);
    c2prn("cap", cap);
    c2prn("car", car);

    Ctxt3P ca3p = mul3(cap, cap, param);
    c3prn("ca3p", ca3p);
    Ctxt3R ca3r = mul3(car, car);
    c3prn("ca3r", ca3r);

    RnsMrs rnsP { 521, 457 };
    EkExtP ekp(skp, param, rsP, rnsP.dynrange_());
    RnsMrs rnsext { rns, rns_ns::Rns::plus, rnsP };
    rns_ns::RnsShrinkRound rshrink(rns, rnsP);
    EkExtR ekr(skr, param, rsR, rnsext, rshrink);

    CtxtP ca2p = relinExt(ca3p, param, ekp);
    CtxtR ca2r = relinExt(ca3r, param, ekr);
    c2prn("ca2p", ca2p);
    c2prn("ca2r", ca2r);

    {
        Poly md2p = decryptP(skp, ca2p, param);
        cout << "md2pU = " << md2p << '\n';
        PolyRns md2r = decryptR(skr, ca2r, param);
        cout << "md2rU = " << md2r << '\n';
        auto a22p = decodeP(param, md2p);
        auto id = param.penc.idelta;
        cout << "a22pU =" << roundv(1e-2, a22p * (1. / id)) << '\n';
        auto a22r = decodeR(param, md2r, rns);
        cout << "a22rU =" << roundv(1e-2, a22r * (1. / id)) << '\n';
    }

    CtxtP ca2scP = rescaleLevel(ca2p, param);
    c2prn("ca2scP", ca2scP);

    Integer qdrop = param.vqs[ca2r.level];
    RnsMrs rnsL {qdrop};
    RnsMrs rnsQ { rns, rns_ns::Rns::minus, rnsL};
    rns_ns::RnsShrinkRound datQ(rnsQ, rnsL);

    CtxtR ca2scR = rescaleLevel(ca2r, datQ);
    c2prn("ca2scR", ca2scR);

    Poly md2p = decryptP(skp, ca2scP, param);
    cout << "md2p = " << md2p << '\n';
    PolyRns md2r = decryptR(skr, ca2scR, param);
    cout << "md2r = " << md2r << '\n';
    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
    auto a22r = decodeR(param, md2r, rns);
    cout << "a22r =" << roundv(1e-2, a22r) << '\n';
}

void t06_mul1()
{
    cout << "\n>>> " << __func__ << " use one function\n";

    using namespace ckks;
    using namespace std::complex_literals;

    ntt::Nttoff nttoff; // FIXME see t05

    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    using rns_ns::RnsMrs;

    Integer delta_(1024);
    Param param(4, Integer(1024), Integer(delta_), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 0.8, 0.5 };

    cout << "a =" << a << '\n';

    RnsMrs rns(param.vqs);
    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    auto c2prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << '\n';
    };

    CtxtP cap = encryptP(skp, map, param, rsP);
    CtxtR car = encryptR(skr, mar, param, rsR);
    c2prn("cap", cap);
    c2prn("car", car);

    RnsMrs rnsP { 521, 457 };
    EkExtP ekp(skp, param, rsP, rnsP.dynrange_());
    RnsMrs rnsext { rns, rns_ns::Rns::plus, rnsP };
    rns_ns::RnsShrinkRound rshrink(rns, rnsP);
    EkExtR ekr(skr, param, rsR, rnsext, rshrink);

    Integer qdrop = param.vqs[car.level];
    RnsMrs rnsL { qdrop };
    RnsMrs rnsQ { rns, rns_ns::Rns::minus, rnsL };
    rns_ns::RnsShrinkRound datQ(rnsQ, rnsL);  // FIXME embed Rns cascade into params

    CtxtP ca2scP = mulExtP(cap, cap, param, ekp);
    CtxtR ca2scR = mulExtR(car, car, param, ekr, datQ);
    c2prn("ca2scP", ca2scP);
    c2prn("ca2scR", ca2scR);

    Poly md2p = decryptP(skp, ca2scP, param);
    cout << "md2p = " << md2p << '\n';
    PolyRns md2r = decryptR(skr, ca2scR, param);
    cout << "md2r = " << md2r << '\n';
    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
    auto a22r = decodeR(param, md2r, rns);
    cout << "a22r =" << roundv(1e-2, a22r) << '\n';
}

void t10_hyb2()
{
    using namespace ckks;
    using namespace std::complex_literals;

    ntt::Nttoff nttoff; // FIXME see t05

    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    using rns_ns::RnsMrs;

    Integer delta_(1024);
    Param param(4, Integer(1024), Integer(delta_), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 0.8, 0.5 };

    cout << "a =" << a << '\n';

    RnsMrs rns(param.vqs);
    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    auto c2prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << '\n';
    };

    CtxtP cap = encryptP(skp, map, param, rsP);
    CtxtR car = encryptR(skr, mar, param, rsR);
    c2prn("cap", cap);
    c2prn("car", car);

    param.w = 16;
    EkHybP ekp(1, skp, param, rsP);

    Integer extP = EkHybR::findExtDigit(rns.getQs(), car.c0.polysize());
    RnsMrs rnsP { extP };
    RnsMrs rnsext { rns, rns_ns::Rns::plus, rnsP };
    rns_ns::RnsShrinkRound rshrink(rns, rnsP);
    ///EkHybR ekr(1, skr, param, rsR);
    //EkHybR ekr(skr, param, rsR, rnsext, rshrink);

    Integer qdrop = param.vqs[car.level];
    RnsMrs rnsL { qdrop };
    RnsMrs rnsQ { rns, rns_ns::Rns::minus, rnsL };
    rns_ns::RnsShrinkRound datQ(rnsQ, rnsL);  // FIXME embed Rns cascade into params

    CtxtP ca2p = mulHybP(cap, cap, param, ekp);
    c2prn("ca2p", ca2p);

    Poly md2p = decryptP(skp, ca2p, param);
    cout << "md2p = " << md2p << '\n';
    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
}
