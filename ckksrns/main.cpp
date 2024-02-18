#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>

#include "ckkselem.h"
#include "ntt.h"
#include "chron.h"
#include "ckkshyb.h"
#include "main.h"

using std::string;
using std::cout;
using std::vector;


using poly::Poly;
using poly::PolyRns;

int main()
try
{
    if (1)
    {
        //t01_rns1();
        //t02_rns2();
        //t03_rns3();
        //t04_rns4();
        /////t04_rnsYes1();
        //t05_rebase();
        //t00_ntt();
        //t01_encode();
        //t02_encSk();
        //t03_encPk();
        //t04_mul3();
        //t04_mul3_v1();
        //t05_mul2_v1();
        //t05_mul2_v2();
        t05_mul2();
        //    t06_mul1();
        //    t07_mul();
    }

    //if (1)
    //{
    //    t08_decomp();
    //    t09_hyb1();
    //    t10_hyb2();
    //}
}

catch (int e)
{
    cout << "error (int): " << e << "\n";
    return 1;
}
catch (string e)
{
    cout << "error (str): " << e << "\n";
    return 1;
}
catch (const char * e)
{
    cout << "error (cc): " << e << "\n";
    return 1;
}
catch (std::exception & e)
{
    cout << "error (std): " << e.what() << "\n";
    return 1;
}
catch (...)
{
    cout << "Unknown exception\n";
    return 1;
}


void t01_rns1()
{
    cout << "\n>>> " << __func__ << '\n';
    using rns_ns::RnsForm;

    rns_ns::RnsMrs rns({13, 5, 7});
    cout << rns.print(); // << '\n';

    int a = 2;
    RnsForm ra(rns, a);
    cout << "ra =" << ra.values() << '\n';
    auto b = ra.lowval();
    cout << "a = " << a << "\tb = " << b << '\n';
}


void t02_rns2()
{
    cout << "\n>>> " << __func__ << '\n';
    using rns_ns::RnsForm;

    rns_ns::RnsMrs rns({3, 5, 7});
    cout << rns.print(); // << '\n';

    int a = 1;
    RnsForm ra1(rns, a);
    cout << "ra =" << ra1.values() << '\n';
    auto ra2 = ra1 + ra1;
    cout << "[+] " << ra2.values() << '\n';
    auto ra4 = ra2 * ra2;
    cout << "[*] " << ra4.values() << '\n';
    auto ra3 = ra4 + - ra1;
    cout << "[~] " << ra3.values() << '\n';
    auto ra = ra3 - ra1;
    cout << "[-] " << ra.values() << '\n';
    auto b = ra.lowval();
    cout << "a = " << a << "\tb = " << b << '\n';
}


void t03_rns3()
{
    cout << "\n>>> " << __func__ << '\n';
    using rns_ns::RnsForm;

    rns_ns::RnsMrs rns({ 11, 5, 7 });
    cout << rns.print(); // << '\n';

    int a = 199;
    RnsForm ra(rns, a);
    cout << "ra =" << ra.values() << '\t' << "(" << ra.blend_() << ")" << '\n';
    while (!ra.islowval())
    {
        ra = ra.pow2div(1);
        cout << "[/2] " << ra.values() << '\t' << "(" << ra.blend_() << ")" << '\n';
    }
    auto b = ra.lowval();
    cout << "a = " << a << "\tb = " << b << '\n';
}

void t04_rns4()
{
    cout << "\n>>> " << __func__ << '\n';
    using rns_ns::RnsForm;

    rns_ns::RnsMrs rns { 5, 7 };
    Integer dyn = rns.dynrange_();
    for (int i = 0; i < dyn; i++)
    {
        for (int j = 1; j < dyn; j++)
        {
            rns_ns::RnsForm a(rns, i), b(rns, j);
            auto fq = a / b;
            auto fr = a % b;
            Integer ifq = fq.blend_();
            Integer ifr = fr.blend_();
            Integer iq = i / j;
            Integer ir = i % j;
            if (ifq != iq) never;
            if (ifr != ir) never;
            //cout << i << '\t' << j << '\t'<< ir << '\t' <<ifr << '\t' << iq << '\t' << ifq << '\n';
        }
    }
}

void t04_rnsYes1()
{
    cout << "\n>>> " << __func__ << '\n';
    using rns_ns::RnsForm;

    //rns_ns::RnsYes rns({ 3, 5, 7 });
    //rns_ns::RnsYes rns({ 17, 23, 29 });
    //rns_ns::RnsYes rns({ 23, 29 });
    //rns_ns::RnsYes rns({ 5, 7 }); // ATTENTION must be 3*Yn << M
    //rns_ns::RnsYes rns({ 7, 9 });
    //rns_ns::RnsYes rns({ 37, 27 });
    //rns_ns::RnsYes rns({ 13, 15, 17 });
    //rns_ns::RnsYes rns({ 13, 23, 17, 19 });
    rns_ns::RnsYes rns({ 9, 5, 7 });
    cout << rns.print(); // << '\n';

    int a = 1 + 0 * 19977;
    RnsForm ra(rns, a);
    cout << "ra =" << ra.values() << '\t' << "(" << ra.blend_() << ")" << '\n';
    while (!ra.islowval())
    {
        ra = ra.pow2div(1);
        cout << "[/2] " << ra.values() << '\t' << "(" << ra.blend_() << ")" << '\n';
    }
    auto b = ra.lowval();
    cout << "a = " << a << "\tb = " << b << '\n';

    Integer dyn = rns.dynrange_();

    if (1)
    {
        RnsForm y = rns.getY();


        for (Integer i = 0; i < dyn; i++)
        {
            RnsForm x(rns, i);
            auto xy = x * y;
            //cout << i << " rho(x)=" << x.rank_() << " rho(x*y)=" << xy.rank_() << '\n';
            cout << i << " " << x.rank_() << " " << xy.rank_() << '\n';
        }
    }

    if (0)
    {
        for (Integer i = 0; i < dyn; i++)
        {
            RnsForm fa(rns, i);
            vint va = fa.values();
            RnsForm fi(rns, rns.negate(va));
            vint vi = fi.values();

            int farho = fa.rank_();
            int firho = fi.rank_();

            auto [xarnk_lo, xarnk_hi, xaapx] = rns.rank_basic_alg(va);
            auto [xirnk_lo, xirnk_hi, xiapx] = rns.rank_basic_alg(vi);

            bool islowa = fa.islowval();
            bool islowi = fi.islowval();
            bool isdiaga = (xarnk_lo != xarnk_hi);
            bool isdiagi = (xirnk_lo != xirnk_hi);

            cout << i << " |" << fa.values() << " | "
                 << "(" << farho << "," << firho << ") "
                 << (islowa ? "L" : "*") << (isdiaga ? "D" : "*")
                 << "[" << xarnk_lo << "," << xarnk_hi << "," << xaapx << "] "
                 << (islowi ? "L" : "*") << (isdiagi ? "D" : "*")
                 << "[" << xirnk_lo << "," << xirnk_hi << "," << xiapx << "]"
                 << '\n';
        }
    }

    if (0)
    {
        for (Integer i = 0 * 3109 + 0 * 27 + 0 * 10505; i < dyn; i++)
        {
            RnsForm ra(rns, i);

            int ri = ra.rank_();
            int pi = i % 2;
            cout << i << '\t' << ri << '\t' << pi << std::flush;

            Integer iblend = ra.blend_();
            if (iblend != i) never;

            //ra.pow2div(1);
            //auto rp = rns.rank_v1(ra.values());
            //auto rp = rns.rank_v2(ra.values());
            //auto rp = rns.rank_v3(ra.values()); // this computes round ksi
            //auto rp4 = rns.rank_v4(ra.values()); // this computes floor ksi - number of fault is more for {13,15,17}
            auto rp3 = rns.rank_v3(ra.values());
            auto rp = rns.rank_v10(ra.values());
            if (rp != rp3) never;
            if (rp.first != ri || rp.second != pi)
            {
                cout << "\tError: computed i,rp: " << i << '\t' << rp.first << '\t' << rp.second;
                ///never;
            }
            cout << '\n';
        }
    }

    if (0)
    {
        for (Integer i = 0 * 27 + 0 * 10505; i < 0 * dyn; i++)
        {
            RnsForm ra(rns, i);
            auto i2 = i / 2;
            cout << i << '\t' << i2 << '\n';
            auto r2 = ra.pow2div(1);
            Integer i3 = r2.blend_();
            //cout << i << '\t' << i2 << '\t' << i3 << '\n';
            if (i2 != i3)
            {
                cout << "Error: i,i2,i3: " << i << '\t' << i2 << '\t' << i3 << '\n';
                never;
            }
        }
    }
}


void t00_ntt()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;

    {
        Poly a;
        a += Integer(1);
        a += Integer(3);
        auto b = a;
        auto c = mul(a, b, Integer(65));
        cout << "a,b,c: " << a << ", " << b << ", " << c << '\n';
    }

    {
        using namespace rns_ns;
        RnsMrs rns {5, 13};
        PolyRns a(rns);
        a += Integer(1);
        a += RnsForm(rns, 3);
        auto b = a;
        auto c = mul(a, b);
        cout << "a,b,c: " << a << ", " << b << ", " << c << '\n';
        cout << "a,b,c: " << a.getCollapsed() << ", " << b << ", " << c << '\n';
    }
}


void t01_encode()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;

    {
        Param param(4, Integer(257), Integer(64), 0);
        cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

        vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

        cout << "a1 =" << a1 << '\n';

        auto ma = encodeP(param, a1);

        cout << "ma = " << ma << '\n';

        auto a2 = decodeP(param, ma);

        cout << "a2 =" << roundv(1e-2, a2) << '\n';
    }

    {
        Param param(4, Integer(257), Integer(64), 0);
        cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

        vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

        cout << "a1 =" << a1 << '\n';

        rns_ns::RnsMrs rns { 15, 17 };

        auto ma = encodeR(param, a1, rns);

        cout << "ma = " << ma << '\n';

        auto a2 = decodeR(param, ma, rns);

        cout << "a2 =" << roundv(1e-2, a2) << '\n';
    }
}

void t02_encSk()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta { 64 };
    ///ParamEncode paramEnc(4, delta);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << param.print() << '\n';
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

    cout << "a1 =" << a1 << '\n';

    rns_ns::RnsMrs rns { 1033, 1009 };
    Poly map = encodeP(param, a1);
    PolyRns mar = encodeR(param, a1, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << map << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    CtxtP cap = encryptP(skp, map, param, rsP);
    cout << "cap.c0 = " << cap.c0 << '\n';
    cout << "cap.c1 = " << cap.c1 << '\n';

    CtxtR car = encryptR(skr, mar, param, rsR);
    cout << "car.c0 = " << car.c0 << '\n';
    cout << "car.c1 = " << car.c1 << '\n';

    Poly m2p = decryptP(skp, cap, param);
    cout << "m2p = " << m2p << '\n';

    PolyRns m2r = decryptR(skr, car, param);
    cout << "m2r = " << m2r << '\n';

    auto a2p = decodeP(param, m2p);
    cout << "a2p =" << roundv(1e-2, a2p) << '\n';
    auto a2r = decodeR(param, m2r, rns);
    cout << "a2r =" << roundv(1e-2, a2r) << '\n';
}



void t03_encPk()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

    cout << "a1 =" << a1 << '\n';

    rns_ns::RnsMrs rns { 1033, 1009 };
    Poly map = encodeP(param, a1);
    PolyRns mar = encodeR(param, a1, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';

    //RndStream rs;
    //Sk sk(rs, param.penc.n);
    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    PkP pkp(skp, param, rsP);
    PkR pkr(skr, param, rsR);

    CtxtP cap = encryptP(pkp, map, param, rsP);
    cout << "cap.c0 = " << cap.c0 << '\n';
    cout << "cap.c1 = " << cap.c1 << '\n';

    CtxtR car = encryptR(pkr, mar, param, rsR);
    cout << "car.c0 = " << car.c0 << '\n';
    cout << "car.c1 = " << car.c1 << '\n';

    poly::Poly m2p = decryptP(skp, cap, param);
    cout << "m2p = " << m2p << '\n';
    poly::PolyRns m2r = decryptR(skr, car, param);
    cout << "m2r = " << m2r << '\n';

    auto a2p = decodeP(param, m2p);
    cout << "a2p =" << roundv(1e-2, a2p) << '\n';
    auto a2r = decodeR(param, m2r, rns);
    cout << "a2r =" << roundv(1e-2, a2r) << '\n';
}

void t04_mul3_v1()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);

    //Param param(2, Integer(1024), Integer(delta), 1);
    //rns_ns::RnsMrs rns{ 1033, 1021 };
    //vector<cx> a = { 3.0 };

    Param param(4, Integer(1024), Integer(delta), 1);
    rns_ns::RnsMrs rns { 1033, 1009 };
    vector<cx> a = { 3.0, 2.0 };

    cout << param.print() << '\n';
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 3.0 + 0.1i, 2.0 + 0i };

    cout << "a =" << a << '\n';

    ///rns_ns::RnsMrs rns { 1033, 1009 };

    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';
    cout << "map decoded =" << roundv(1e-2, decodeP(param, map)) << '\n';
    cout << "mar decoded =" << roundv(1e-2, decodeR(param, mar, rns)) << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    CtxtP cap = encryptP(skp, map, param, rsP);
    cout << "cap.c0 = " << cap.c0 << '\n';
    cout << "cap.c1 = " << cap.c1 << '\n';

    CtxtR car = encryptR(skr, mar, param, rsR);
    cout << "car.c0 = " << car.c0 << '\n';
    cout << "car.c1 = " << car.c1 << '\n';

    {
        Poly mdp = decryptP(skp, cap, param);
        cout << "mdp = " << mdp << '\n';
        auto a2p = decodeP(param, mdp);
        cout << "a2p =" << roundv(1e-2, a2p) << '\n';

        PolyRns mdr = decryptR(skr, car, param);
        cout << "mdr = " << mdr << '\n';
        auto a2r = decodeR(param, mdr, rns);
        cout << "a2r =" << roundv(1e-2, a2r) << '\n';
    }

    Ctxt3P ca3p = mul3(cap, cap, param);
    Ctxt3R ca3r = mul3(car, car);

    auto c3prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << c.c2 << '\n';
    };

    c3prn("ca3p", ca3p);
    c3prn("ca3r", ca3r);

    {
        Poly md3p = decryptP3(skp, ca3p, param);
        cout << "md3p = " << md3p << '\n';
        PolyRns md3r = decryptR3(skr, ca3r, param);
        cout << "md3r = " << md3r << '\n';

        auto a23p = decodeP(param, md3p);
        cout << "a23p =" << roundv(1e-2, a23p) << '\n';
        auto a23r = decodeR(param, md3r, rns);
        cout << "a23r =" << roundv(1e-2, a23r) << '\n';

        Poly md3scP = rescaleRound(md3p, param.penc.idelta);
        auto a23scP = decodeP(param, md3scP);
        cout << "a23scP =" << roundv(1e-4, a23scP) << '\n';

        PolyRns md3scR = rescaleRoundRns(md3r, param.penc.idelta);
        auto a23scR = decodeR(param, md3scR, rns);
        cout << "a23scR =" << roundv(1e-4, a23scR) << '\n';
    }

    cout << '\n';
    {
        Ctxt3P ca3scP = rescale(ca3p, param.penc.idelta, param);
        c3prn("ca3scP", ca3scP);
        Ctxt3R ca3scR = rescale(ca3r, param.penc.idelta);
        c3prn("ca3scR", ca3scR);

        Poly md3sp = decryptP3(skp, ca3scP, param);
        cout << "md3sp = " << md3sp << '\n';
        PolyRns md3sr = decryptR3(skr, ca3scR, param);
        cout << "md3sr = " << md3sr << '\n';

        auto a23sp = decodeP(param, md3sp);
        cout << "a23sp =" << roundv(1e-4, a23sp) << '\n';
        auto a23sr = decodeR(param, md3sr, rns);
        cout << "a23sr =" << roundv(1e-4, a23sr) << '\n';
    }
}

void t04_mul3()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);
    //Param param(2, Integer(1024), Integer(delta), 1);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << param.print() << '\n';
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    //vector<cx> a = { 3.0 };
    vector<cx> a = { 3.0, 2.0 };
    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 3.0 + 0.1i, 2.0 + 0i };

    cout << "a =" << a << '\n';

    rns_ns::RnsMrs rns { 1033, 1009 };

    Poly map = encodeP(param, a);
    PolyRns mar = encodeR(param, a, rns);

    cout << "map = " << map << '\n';
    cout << "mar = " << mar << '\n';
    cout << "map decoded =" << roundv(1e-2, decodeP(param, map)) << '\n';
    cout << "mar decoded =" << roundv(1e-2, decodeR(param, mar, rns)) << '\n';

    RndStream rsP, rsR;
    SkP skp(rsP, param.penc.n);
    SkR skr(rsR, param.penc.n, rns);

    CtxtP cap = encryptP(skp, map, param, rsP);
    cout << "cap.c0 = " << cap.c0 << '\n';
    cout << "cap.c1 = " << cap.c1 << '\n';

    CtxtR car = encryptR(skr, mar, param, rsR);
    cout << "car.c0 = " << car.c0 << '\n';
    cout << "car.c1 = " << car.c1 << '\n';

    {
        Poly mdp = decryptP(skp, cap, param);
        cout << "mdp = " << mdp << '\n';
        auto a2p = decodeP(param, mdp);
        cout << "a2p =" << roundv(1e-2, a2p) << '\n';

        PolyRns mdr = decryptR(skr, car, param);
        cout << "mdr = " << mdr << '\n';
        auto a2r = decodeR(param, mdr, rns);
        cout << "a2r =" << roundv(1e-2, a2r) << '\n';
    }

    {
        Poly ma2p = add(map, map, param.q0());
        cout << "xa2addp =" << roundv(1e-2, decodeP(param, ma2p)) << '\n';

        PolyRns ma2r = add(mar, mar);
        cout << "xa2addr =" << roundv(1e-2, decodeR(param, ma2r, rns)) << '\n';
    }

    if (0)
    {
        cout << "\nWierd calc - attempt to match, but it doesnt make sense\n";
        Poly ma2p = mul_simple(map, map);
        cout << "ma2p = " << ma2p << '\n';
        cout << "xa2mulp =" << roundv(1e-2, decodeP(param, ma2p)) << '\n';
        Poly ma2scp = rescaleRound(ma2p, param.penc.idelta);
        cout << "ma2mscp =" << roundv(1e-2, decodeP(param, ma2scp)) << '\n';
        PolyRns ma2r = mul_simple(mar, mar);
        cout << "ma2r = " << ma2r << '\n';
        cout << "xa2mulr =" << roundv(1e-2, decodeR(param, ma2r, rns)) << '\n';
    }

    //never;

    cout << "\nadd\n";
    {
        CtxtP ca2p = add(cap, cap, param);
        Poly md2p = decryptP(skp, ca2p, param);
        cout << "md2p = " << md2p << '\n';
        auto a22p = decodeP(param, md2p);
        cout << "a22p =" << roundv(1e-2, a22p) << '\n';

        CtxtR ca2r = add(car, car);
        PolyRns md2r = decryptR(skr, ca2r, param);
        cout << "md2r = " << md2r << '\n';
        auto a22r = decodeR(param, md2r, rns);
        cout << "a22 =" << roundv(1e-2, a22r) << '\n';
    }

    //never;

    cout << "\nmul" << '\n';
    Ctxt3P ca3p = mul3(cap, cap, param);
    Ctxt3R ca3r = mul3(car, car);

    auto c3prn = [](string nm, auto c)
    {
        cout << nm << " = " << c.c0 << c.c1 << c.c2 << '\n';
    };

    c3prn("ca3p", ca3p);
    c3prn("ca3r", ca3r);

    {
        Poly md3p = decryptP3(skp, ca3p, param);
        cout << "md3p = " << md3p << '\n';
        PolyRns md3r = decryptR3(skr, ca3r, param);
        cout << "md3r = " << md3r << '\n';

        auto a23p = decodeP(param, md3p);
        cout << "a23p =" << roundv(1e-2, a23p) << '\n';
        auto a23r = decodeR(param, md3r, rns);
        cout << "a23r =" << roundv(1e-2, a23r) << '\n';

        Poly md3scP = rescaleRound(md3p, param.penc.idelta);
        auto a23scP = decodeP(param, md3scP);
        cout << "a23scP =" << roundv(1e-2, a23scP) << '\n';

        PolyRns md3scR = rescaleRoundRns(md3r, param.penc.idelta);
        auto a23scR = decodeR(param, md3scR, rns);
        cout << "a23scR =" << roundv(1e-2, a23scR) << '\n';
    }

    cout << '\n';
    {
        Ctxt3P ca3scP = rescale(ca3p, param.penc.idelta, param);
        c3prn("ca3scP", ca3scP);
        Ctxt3R ca3scR = rescale(ca3r, param.penc.idelta);
        c3prn("ca3scR", ca3scR);

        Poly md3sp = decryptP3(skp, ca3scP, param);
        cout << "md3sp = " << md3sp << '\n';
        PolyRns md3sr = decryptR3(skr, ca3scR, param);
        cout << "md3sr = " << md3sr << '\n';

        auto a23sp = decodeP(param, md3sp);
        cout << "a23sp =" << roundv(1e-2, a23sp) << '\n';
        auto a23sr = decodeR(param, md3sr, rns);
        cout << "a23sr =" << roundv(1e-2, a23sr) << '\n';
    }
}

