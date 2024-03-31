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




void t10_hyb2_b1()
{
    using namespace ckks;
    using namespace std::complex_literals;

    ntt::Nttoff nttoff; // FIXME see t05

    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    using rns_ns::RnsMrs;

    Integer delta_(1024);
    Param param(2, Integer(1024), Integer(delta_), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 0.80 };
    //vector<cx> a = { 0.8, 0.5 };

    cout << "a =" << a << '\n';

    RnsMrs rns(param.vqs);
    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';

    RndStream rsP, rsR, rsE;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    auto c2prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << '\n';
    };

    CtxtP cap = encryptP(skp, map, param, rsP);
    CtxtR car = encryptR(skr, mar, param, rsR);
    CtxtR cae = encryptR(skr, mar, param, rsE); // just to use rsE
    c2prn("cap", cap);
    c2prn("car", car);

    param.w = 16;
    EkHybP ekp(1, skp, param, rsP);

    Integer extP = EkHybR::findExtDigit(rns.getQs(), car.c0.polysize());
    RnsMrs rnsP { extP };
    RnsMrs rnsext { rns, rns_ns::Rns::plus, rnsP };
    rns_ns::RnsShrinkRound rshrink(rns, rnsP);
    ///EkHybR ekr(1, skr, param, rsR);
    EkHybR ekr(skr, param, rsR, rnsext, rshrink);
    EkExtR eke(skr, param, rsE, rnsext, rshrink);

    Integer qdrop = param.vqs[car.level];
    RnsMrs rnsL { qdrop };
    RnsMrs rnsQ { rns, rns_ns::Rns::minus, rnsL };
    rns_ns::RnsShrinkRound datQ(rnsQ, rnsL);  // FIXME embed Rns cascade into params

    CtxtP ca2p = mulHybP(cap, cap, param, ekp);
    c2prn("ca2p", ca2p);

    CtxtR ca2r = mulHybR(car, car, param, ekr, datQ);
    c2prn("ca2r", ca2r);

    CtxtR ca2e = mulExtR(car, car, param, eke, datQ);
    c2prn("ca2e", ca2e);

    Poly md2p = decryptP(skp, ca2p, param);
    //md2p = rangeUpP(md2p, 1024); // FIXME cludge
    cout << "md2p = " << md2p << " (1-x)"  << '\n';
    auto md2r = decryptR(skr, ca2r, param);
    cout << "md2r = " << md2r << '\n';
    auto md2e = decryptR(skr, ca2e, param);
    cout << "md2e = " << md2e << '\n';

    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
    auto a22r = decodeR(param, md2r, rns);
    cout << "a22r =" << roundv(1e-2, a22r) << '\n';
    auto a22e = decodeR(param, md2e, rns);
    cout << "a22e =" << roundv(1e-2, a22e) << '\n';
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
    EkHybR ekr(skr, param, rsR, rnsext, rshrink);

    Integer qdrop = param.vqs[car.level];
    RnsMrs rnsL { qdrop };
    RnsMrs rnsQ { rns, rns_ns::Rns::minus, rnsL };
    rns_ns::RnsShrinkRound datQ(rnsQ, rnsL);  // FIXME embed Rns cascade into params

    CtxtP ca2p = mulHybP(cap, cap, param, ekp);
    c2prn("ca2p", ca2p);

    CtxtR ca2r = mulHybR(car, car, param, ekr, datQ);
    c2prn("ca2r", ca2r);

    Poly md2p = decryptP(skp, ca2p, param);
    cout << "md2p = " << md2p << '\n';
    auto md2r = decryptR(skr, ca2r, param);
    cout << "md2r = " << md2r << '\n';

    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
    auto a22r = decodeR(param, md2r, rns);
    cout << "a22r =" << roundv(1e-2, a22r) << '\n';

    if (1) // FastBConv
    {
        CtxtR ca2f = mulHybR_fbc(car, car, param, ekr, datQ);
        c2prn("ca2f", ca2f);

        auto md2f = decryptR(skr, ca2f, param);
        cout << "md2f = " << md2f << '\n';

        auto a22f = decodeR(param, md2f, rns);
        cout << "a22f =" << roundv(1e-2, a22f) << '\n';
    }
}

