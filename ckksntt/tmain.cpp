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


void ntt_timing_run(int n, int mod)
{
    Poly b(n, Integer(1));

    auto btime = chron::now();
    for ( int i = 0; i < 1000; i++ )
        b = ntt::nttBfly(b, Integer(mod));
    cout << n << ": " << (chron::now() - btime) / 1000.0 << "ms\n";

    std::ignore = b;
}

void ntt_timing()
{
    /*
    512/17921
    1024/39937
    2048/83969
    4096/151553
    8192/286721
    16384/638977
    32768/1146881
    65536/2424833
    */
    ntt_timing_run(512, 39937);
    ntt_timing_run(1024, 83969);
    ntt_timing_run(2048, 151553);
    ntt_timing_run(4096, 286721);
    ntt_timing_run(8192, 638977);
    ntt_timing_run(16384, 1146881);
    ntt_timing_run(32768, 2424833);
}

void t00_ntt()
{
    using namespace ckks;
    Poly a;
    a += Integer(1);
    a += Integer(3);
    auto b = a;
    auto c = mul(a, b, Integer(5));
    //cout << "a,b,c: " << a << ", " << b << ", " << c << '\n';
}

void t01_encode()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Param param(4, Integer(257), Integer(64), 0);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

    cout << "a1 =" << a1 << '\n';

    auto ma = encode(param, a1);

    cout << "ma = " << ma << '\n';

    auto a2 = decode(param, ma);

    cout << "a2 =" << roundv(1e-2, a2) << '\n';
}

void t02_encSk()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta {64};
    ///ParamEncode paramEnc(4, delta);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

    cout << "a1 =" << a1 << '\n';

    Poly ma = encode(param, a1);

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    ///ParamQ pq(Integer(1024), 1, paramEnc);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    Poly m2 = decrypt(sk, ca, param);
    cout << "m2 = " << m2 << '\n';

    auto a2 = decode(param, m2);

    cout << "a2 =" << roundv(1e-2, a2) << '\n';
}

void t03_encPk()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);
    ///ParamEncode paramEnc(4, delta);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a1 = { 3.0 + 4i, 2.0 - 1i };

    cout << "a1 =" << a1 << '\n';

    Poly ma = encode(param, a1);

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    ///ParamQ pq(Integer(1024), 1, paramEnc);

    Pk pk(sk, param, rs);

    Ctxt ca = encrypt(pk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    poly::Poly m2 = decrypt(sk, ca, param);
    cout << "m2 = " << m2 << '\n';

    auto a2 = decode(param, m2);

    cout << "a2 =" << roundv(1e-2, a2) << '\n';
}

void t04_mul3()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a = { 3.0, 2.0 };
    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 3.0 + 0.1i, 2.0 + 0i };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);

    cout << "ma = " << ma << '\n';
    cout << "ma decoded =" << roundv(1e-2, decode(param, ma)) << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    ///ParamQ pq(Integer(1024), 1, paramEnc);

    //Pk pk(sk, param, rs);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    //Ctxt3 ca3 = mul3(ca,ca);

    {
        Poly md = decrypt(sk, ca, param);
        cout << "md = " << md << '\n';
        auto a2 = decode(param, md);
        cout << "a2 =" << roundv(1e-2, a2) << '\n';
    }

    {
        Poly ma2 = add(ma, ma, param.q0());
        cout << "xa2add =" << roundv(1e-2, decode(param, ma2)) << '\n';
    }

    {
        Poly ma2 = mul_simple(ma, ma);
        Poly ma2sc = rescaleRound(ma2, param.penc.idelta);
        cout << "xa2mul =" << roundv(1e-2, decode(param, ma2sc)) << '\n';
    }

    {
        Ctxt ca2 = add(ca, ca, param);
        Poly md2 = decrypt(sk, ca2, param);
        cout << "md2 = " << md2 << '\n';
        auto a22 = decode(param, md2);
        cout << "a22 =" << roundv(1e-2, a22) << '\n';
    }

    cout << '\n';
    Ctxt3 ca3 = mul3(ca, ca, param);
    {
        Poly md3 = decrypt(sk, ca3, param);
        cout << "md3 = " << md3 << '\n';
        auto a23 = decode(param, md3);
        cout << "a23 =" << roundv(1e-2, a23) << '\n';

        Poly md3sc = rescaleRound(md3, param.penc.idelta);
        auto a23sc = decode(param, md3sc);
        cout << "a23sc =" << roundv(1e-2, a23sc) << '\n';
    }

    cout << '\n';
    {
        Ctxt3 ca3sc = rescale(ca3, param.penc.idelta);
        Poly md3s = decrypt(sk, ca3sc, param);
        cout << "md3s = " << md3s << '\n';
        auto a23s = decode(param, md3s);
        cout << "a23s =" << roundv(1e-2, a23s) << '\n';
    }
}

void t05_mul2()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);

    vector<cx> a = { 3.0, 2.0 };
    //vector<cx> a = { 1.0, 0.0 };
    //vector<cx> a = { 1.5625, 1.5625 };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);

    if (0)
    {
        ma.v[0] = ma.v[1] = ma.v[2] = ma.v[3] = Integer(0);
        ma.v[0] = Integer(100);
        cout << "ma xx = " << ma << '\n';
        auto a2 = decode(param, ma);
        cout << "a2 xx = " << roundv(1e-16, a2) << '\n';
        return;
    }

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    Ctxt3 ca3 = mul3(ca, ca, param);
    cout << "ca3.c0 = " << ca3.c0 << '\n';
    cout << "ca3.c1 = " << ca3.c1 << '\n';
    cout << "ca3.c2 = " << ca3.c2 << '\n';
    EkExt ek(sk, param, rs);
    Ctxt ca2 = relinExt(ca3, param, ek);
    cout << "ca2.c0 = " << ca2.c0 << '\n';
    cout << "ca2.c1 = " << ca2.c1 << '\n';
    Ctxt ca2sc = rescale(ca2, param.penc.idelta);
    cout << "ca2sc.c0 = " << ca2sc.c0 << '\n';
    cout << "ca2sc.c1 = " << ca2sc.c1 << '\n';

    Poly md2 = decrypt(sk, ca2sc, param);
    cout << "md2 = " << md2 << '\n';
    auto a22 = decode(param, md2);
    cout << "a22 =" << roundv(1e-2, a22) << '\n';
}

void t06_mul1()
{
    cout << "\n>>> " << __func__ << " use one function\n";

    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);

    vector<cx> a = { 3.0, 2.0 };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);
    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    Ctxt ca = encrypt(sk, ma, param, rs);
    EkExt ek(sk, param, rs);
    Ctxt ca2sc = mulExt(ca, ca, param, ek);

    Poly md2 = decrypt(sk, ca2sc, param);
    cout << "md2 = " << md2 << '\n';
    auto a22 = decode(param, md2);
    cout << "a22 =" << roundv(1e-2, a22) << '\n';
}

void t07_mul()
{
    cout << "\n>>> " << __func__ << " L=2\n";

    using namespace ckks;
    using namespace std::complex_literals;
}


void t08_decomp()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);
    cout << param.print() << '\n';

    vector<cx> a = { 0.5, 2.0 };
    //vector<cx> b = { 4.0, 3.0 };
    vector<cx> b = { 3.0, 1.0 };
    //vector<cx> a = { 1.5625, 1.5625 };

    cout << "a =" << a << '\n';
    cout << "b =" << b << '\n';

    Poly ma = encode(param, a);
    Poly mb = encode(param, b);

    cout << "ma = " << ma << '\n';
    cout << "mb = " << mb << '\n';

    {
        cout << "\nsimple\n";
        Poly mcSc = poly::mul_simple(ma, mb);
        Poly mc = rescaleRound(mcSc, param.penc.idelta);

        cout << "mc = ma*mb = " << mc << '\n';

        vector<cx> c = decode(param, mc);
        cout << "c =" << roundv(1e-2, c) << '\n';
    }

    auto qL = param.qL();

    {
        cout << "\nrangeUp\n";
        auto maU = rangeUp(ma, qL);
        auto mbU = rangeUp(mb, qL);
        cout << "maU = " << maU << '\n';
        cout << "mbU = " << mbU << '\n';
        Poly mcSc = poly::mul(maU, mbU, qL);
        cout << "mcSc = " << mcSc << '\n';
        Poly mcQ = rescaleRound(mcSc, param.penc.idelta);
        cout << "mcQ = " << mcQ << '\n';
        Poly mc = rangeCenter(mcQ, param.q(param.levels - 1));

        cout << "mc = ma*mb = " << mc << '\n';

        vector<cx> c = decode(param, mc);
        cout << "c =" << roundv(1e-2, c) << '\n';
    }

    param.w = 16;

    {
        cout << "\ndecomposition\n";
        auto maU = rangeUp(ma, qL);
        auto mbU = rangeUp(mb, qL);
        Poly mcSc = poly::mul(maU, mbU, qL);
        cout << "mcSc = " << mcSc << '\n';
        Poly mcQ = rescaleRound(mcSc, param.penc.idelta);
        cout << "mcQ = " << mcQ << '\n';
        Poly mc = rangeCenter(mcQ, param.q(param.levels - 1));
        cout << "mc = ma*mb = " << mc << '\n';
        vector<cx> c = decode(param, mc);
        cout << "c =" << roundv(1e-2, c) << '\n';

        cout << '\n';
        if (param.w == Integer(0)) never;

        auto wda = poly::WD(maU, param.w, qL);
        cout << "wda = \n";
        for (auto x : wda) cout << " " << x << '\n';

        auto pwb = poly::PW(mbU, param.w, qL);
        cout << "pwb = \n";
        for (auto x : pwb) cout << " " << x << '\n';

        Poly ab = poly::dot(wda, pwb, qL);
        cout << "ab = " << ab << '\n';
        cout << "<> = " << mcSc << '\n';
    }

}

void t09_hyb1()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);

    vector<cx> a = { 3.0, 2.0 };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    Ctxt3 ca3 = mul3(ca, ca, param);
    cout << "ca3.c0 = " << ca3.c0 << '\n';
    cout << "ca3.c1 = " << ca3.c1 << '\n';
    cout << "ca3.c2 = " << ca3.c2 << '\n';

    param.w = 16;
    EkHyb ek(1, sk, param, rs);
    cout << "\nek=\n" << ek.print() << '\n';
    Ctxt ca2 = ckks::relinHyb(ca3, param, ek);
    cout << "ca2.c0 = " << ca2.c0 << '\n';
    cout << "ca2.c1 = " << ca2.c1 << '\n';
    Ctxt ca2sc = rescale(ca2, param.penc.idelta);
    cout << "ca2sc.c0 = " << ca2sc.c0 << '\n';
    cout << "ca2sc.c1 = " << ca2sc.c1 << '\n';

    Poly md2 = decrypt(sk, ca2sc, param);
    cout << "md2 = " << md2 << '\n';
    auto a22 = decode(param, md2);
    cout << "a22 =" << roundv(1e-2, a22) << '\n';
}

void t10_hyb2()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param(4, Integer(1024), Integer(delta), 1);

    vector<cx> a = { 3.0, 2.0 };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    param.w = 16;
    EkHyb ek(1, sk, param, rs);

    Ctxt ca2 = mulHyb(ca, ca, param, ek);
    cout << "ca2.c0 = " << ca2.c0 << '\n';
    cout << "ca2.c1 = " << ca2.c1 << '\n';

    Poly md2 = decrypt(sk, ca2, param);
    cout << "md2 = " << md2 << '\n';
    auto a22 = decode(param, md2);
    cout << "a22 =" << roundv(1e-2, a22) << '\n';
}


void t04_mul3_b1()
{
    cout << "\n>>> " << __func__ << '\n';

    using namespace ckks;
    using namespace std::complex_literals;
    Integer delta(64);
    Param param(2, Integer(1024), Integer(delta), 1);
    cout << "xi =" << roundv(1e-3, param.penc.vxi) << '\n';

    vector<cx> a = { 3.0 };

    cout << "a =" << a << '\n';

    Poly ma = encode(param, a);

    cout << "ma = " << ma << '\n';
    cout << "ma decoded =" << roundv(1e-2, decode(param, ma)) << '\n';

    RndStream rs;
    Sk sk(rs, param.penc.n);

    Ctxt ca = encrypt(sk, ma, param, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    {
        Poly md = decrypt(sk, ca, param);
        cout << "md = " << md << '\n';
        auto a2 = decode(param, md);
        cout << "a2 =" << roundv(1e-2, a2) << '\n';
    }

    {
        Poly ma2 = add(ma, ma, param.q0());
        cout << "xa2add =" << roundv(1e-2, decode(param, ma2)) << '\n';
    }

    {
        Poly ma2 = mul_simple(ma, ma);
        Poly ma2sc = rescaleRound(ma2, param.penc.idelta);
        cout << "xa2mul =" << roundv(1e-2, decode(param, ma2sc)) << '\n';
    }

    {
        Ctxt ca2 = add(ca, ca, param);
        Poly md2 = decrypt(sk, ca2, param);
        cout << "md2 = " << md2 << '\n';
        auto a22 = decode(param, md2);
        cout << "a22 =" << roundv(1e-2, a22) << '\n';
    }

    cout << '\n';
    Ctxt3 ca3 = mul3(ca, ca, param);
    {
        Poly md3 = decrypt(sk, ca3, param);
        cout << "md3 = " << md3 << '\n';
        auto a23 = decode(param, md3);
        cout << "a23 =" << roundv(1e-2, a23) << '\n';

        Poly md3sc = rescaleRound(md3, param.penc.idelta);
        auto a23sc = decode(param, md3sc);
        cout << "a23sc =" << roundv(1e-2, a23sc) << '\n';
    }

    cout << '\n';
    {
        Ctxt3 ca3sc = rescale(ca3, param.penc.idelta);
        Poly md3s = decrypt(sk, ca3sc, param);
        cout << "md3s = " << md3s << '\n';
        auto a23s = decode(param, md3s);
        cout << "a23s =" << roundv(1e-2, a23s) << '\n';
    }
}


void t05_mul2_b1()
{
    cout << "\n>>> " << __func__ << '\n';

    struct Nttoff
    {
        bool oldval;
        Nttoff() : oldval(ntt::disabled) { ntt::disabled = true; }
        ~Nttoff() { ntt::disabled = oldval; }
    } nttoff;


    using namespace ckks;
    using namespace std::complex_literals;

    Integer delta(64);
    Param param2(2, Integer(1024), Integer(delta), 1);
    Param param1(1, Integer(1024), Integer(delta), 1);

    vector<cx> a = { 3.0 };

    cout << "a =" << a << '\n';

    Poly ma0 = encode(param2, a);
    Poly ma; ma += ma0.v[0];

    if (0)
    {
        ma.v[0] = ma.v[1] = ma.v[2] = ma.v[3] = Integer(0);
        ma.v[0] = Integer(100);
        cout << "ma xx = " << ma << '\n';
        auto a2 = decode(param2, ma);
        cout << "a2 xx = " << roundv(1e-16, a2) << '\n';
        return;
    }

    cout << "ma = " << ma << '\n';

    RndStream rs;
    Sk sk(rs, param1.penc.n);

    Ctxt ca = encrypt(sk, ma, param1, rs);
    cout << "ca.c0 = " << ca.c0 << '\n';
    cout << "ca.c1 = " << ca.c1 << '\n';

    {
        Poly md = decrypt(sk, ca, param1);
        cout << "md = " << md << '\n';
        Poly md0 = md; md0 += Integer(0);
        auto a2 = decode(param2, md0);
        cout << "a2 =" << roundv(1e-2, a2) << '\n';
    }

    Ctxt3 ca3 = mul3(ca, ca, param1);

    {
        Poly md3 = decrypt(sk, ca3, param1);
        cout << "md3 = " << md3 << '\n';
        Poly md30 = md3; md30 += Integer(0);
        auto a23 = decode(param2, md30);
        cout << "a23 =" << roundv(1e-2, a23) << '\n';

        Poly md3sc = rescaleRound(md3, param1.penc.idelta);
        Poly md3sc0 = md3sc; md3sc0 += Integer(0);
        auto a23sc = decode(param2, md3sc0);
        cout << "a23sc =" << roundv(1e-2, a23sc) << '\n';
    }

    cout << "ca3.c0 = " << ca3.c0 << '\n';
    cout << "ca3.c1 = " << ca3.c1 << '\n';
    cout << "ca3.c2 = " << ca3.c2 << '\n';
    EkExt ek(sk, param1, rs);
    Ctxt ca2 = relinExt(ca3, param1, ek);
    cout << "ca2.c0 = " << ca2.c0 << '\n';
    cout << "ca2.c1 = " << ca2.c1 << '\n';
    Ctxt ca2sc = rescale(ca2, param1.penc.idelta);
    cout << "ca2sc.c0 = " << ca2sc.c0 << '\n';
    cout << "ca2sc.c1 = " << ca2sc.c1 << '\n';

    Poly md2 = decrypt(sk, ca2sc, param1);
    cout << "md2 = " << md2 << '\n';
    Poly md20 = md2; md20 += Integer(0);
    auto a22 = decode(param2, md20);
    cout << "a22 =" << roundv(1e-2, a22) << '\n';
}

