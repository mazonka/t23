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
    struct Nttoff
    {
        bool oldval;
        Nttoff() : oldval(ntt::disabled) { ntt::disabled = true; }
        ~Nttoff() { ntt::disabled = oldval; }
    } nttoff;

    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    using rns_ns::RnsMrs;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 3.0, 2.0 };
    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 1.5625, 1.5625 };

    cout << "a =" << a << '\n';

    RnsMrs rns { 1033, 1009 };
    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';

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

    CtxtP ca2scP = rescale(ca2p, param.penc.idelta, param);
    c2prn("ca2scP", ca2scP);
    //CtxtR ca2scR = rescale(ca2r, param.penc.idelta, param);
    //c2prn("ca2scR", ca2scR);

    Poly md2p = decryptP(skp, ca2scP, param);
    cout << "md2p = " << md2p << '\n';
    //PolyRns md2r = decryptR(skr, ca2scR, param);
    //cout << "md2r = " << md2r << '\n';
    auto a22p = decodeP(param, md2p);
    cout << "a22p =" << roundv(1e-2, a22p) << '\n';
    //auto a22r = decodeR(param, md2r, rns);
    //cout << "a22r =" << roundv(1e-2, a22r) << '\n';
}

void t05_rebase()
{
    int step = 4;

    using namespace rns_ns;
    RnsMrs rnsa { 3, 5, 7 };
    RnsMrs rnsb { 7, 9, 11 };

    auto dyna = rnsa.dynrange_();

    for (Integer i = 0; i < dyna; i += dyna / step)
    {
        RnsForm x(rnsa, i);
        auto y = x.rebaseAny(rnsb);
        auto z = y.rebaseAny(rnsa);
        if (z != x) never;
        Integer iz = z.blend_();
        if (iz != i) never;
        cout << i << '\t' << x.values() << '\t' << y.values() << '\t' << z.values() << '\n';
    }

    RnsMrs rnsc { 7, 9, 11, 5, 13, 17 };
    auto dynb = rnsb.dynrange_();

    for (Integer i = 0; i < dynb; i += dynb / step)
    {
        RnsForm x(rnsb, i);
        auto y = x.rebaseAny(rnsc);
        auto z = y.rebaseCut(rnsb);
        if (z != x) never;
        Integer iz = z.blend_();
        if (iz != i) never;
        cout << i << '\t' << x.values() << '\t' << y.values() << '\t' << z.values() << '\n';
    }

    {
        RnsMrs rnsQ { 7 };
        RnsMrs rnsQP = rnsa;
        RnsMrs rnsP(rnsQP, Rns::minus, rnsQ);

        double ddynP = rnsP.dynrangeDbl();
        RnsShrinkRound data1(rnsQ, rnsP);
        for (Integer i = 0; i < dyna; i += 31)
        {
            RnsForm x(rnsQP, i);
            auto y = x.rebaseShrinkRound(data1);
            Integer iy = y.blend_();
            //if (iz != i) never;
            cout << i << '\t' << iy << '\t' << (i / ddynP) << '\t'
                 << x.values() << '\t' << y.values() << '\n';
        }

        double ddynQ = rnsQ.dynrangeDbl();
        RnsShrinkRound data2(rnsP, rnsQ);
        for (Integer i = 0; i < dyna; i += 31)
        {
            RnsForm x(rnsQP, i);
            auto y = x.rebaseShrinkRound(data2);
            Integer iy = y.blend_();
            //if (iz != i) never;
            cout << i << '\t' << iy << '\t' << (i / ddynQ) << '\t'
                 << x.values() << '\t' << y.values() << '\n';
        }
    }

    {
        RnsMrs rnsP { 7, 9, 11 };
        RnsMrs rnsQ { 5, 17, 13 };
        RnsMrs rnsQP(rnsP, Rns::plus, rnsQ);

        double ddynP = rnsP.dynrangeDbl();
        Integer dynQ = rnsQ.dynrange_();
        Integer dynQP = rnsQP.dynrange_();
        RnsShrinkRound data1(rnsQ, rnsP);
        for (Integer i = 0; i < dynQP; i += 100000)
        {
            RnsForm x(rnsQP, i);
            auto y = x.rebaseShrinkRound(data1);
            Integer iy = y.blend_();
            double r = (i / ddynP);

            cout << i << '\t' << iy << '\t' << x.values() << '\t'
                 << y.values() << '\t' << r << '\n';

            auto z = (Integer)(r + 0.5);
            auto w = z % dynQ;
            if (iy != w ) never;
        }
    }

    {
        auto dync = rnsc.dynrange_();

        //double cbP = 1.0 * dync / dynb;
        double cbP = 1.0 * dynb;

        RnsMrs rnsCmB(rnsc, Rns::minus, rnsb);
        RnsShrinkRound datc(rnsCmB, rnsb);

        for (Integer i = 0; i < dync; i += 100027)
        {
            RnsForm x(rnsc, i);
            auto y = x.rebaseShrinkRound(datc);
            Integer iy = y.blend_();
            //if (iz != i) never;
            cout << i << '\t' << iy << '\t' << x.values() << '\t'
                 << y.values() << '\t' << (i / cbP) << '\n';
        }
    }
}

