// ref: Mikhail Selianinau, 2021 Computationally Efficient Approach ...
#include <algorithm>

#include "ccrun.h"
//#include "ccrut.h"

using Integer = int;
using RnsForm = std::vector<Integer>;







template <class T>
T egcdT_rec(const T & a, const T & b, T & s0, T & t0, T s1, T t1, T s2, T t2)
{
    auto c = a % b;
    if (c == T(0)) return b;
    auto q = (a - c) / b; // exact division
    s0 = s2 - s1 * q;
    t0 = t2 - t1 * q;
    return egcdT_rec(b, c, s0, t0, s0, t0, s1, t1);
}

template <class T>
T egcdT(const T & a, const T & b, T & s0, T & t0)
{
    T s1(0), t1(1), s2(1), t2(0);
    return egcdT_rec(a, b, s0, t0, s1, t1, s2, t2);
}

Integer invert(Integer a, Integer mod)
{
    Integer s, t;
    auto gcd = egcdT<Integer>(a, mod, s, t);
    if (gcd * gcd != 1 ) nevers("no inverse (ull)");
    while (s < 0) s += mod;
    s /= gcd;
    return s % mod;
}



using vint = std::vector<Integer>;

struct Rns
{
        std::vector<Integer> qs, Ms, us, Bs;
        using vvint = std::vector<vint>;
        vvint cjk;
        Integer Y;


        RnsForm split(Integer x) const;
        Integer blend(const RnsForm & v) const;

        void init(std::initializer_list<Integer> li);
        Rns(std::initializer_list<Integer> li) { init(li); }

        Integer dynrange() const { Integer r = 1; for (auto x : qs) r *= x; return r; }
        string print() const;
        int size() const { return (int)qs.size(); }
        Integer EBx(RnsForm x) const;
        Integer EMx(RnsForm x) const;
        vint mrs(RnsForm x) const;
        vint xi(RnsForm x) const;

    private:
        void assertRnsInited() const { if (qs.empty()) nevers("rns is not set"); }
};

void Rns::init(std::initializer_list<Integer> li)
{
    qs = li;
    std::sort(qs.begin(), qs.end());
    Integer d = dynrange();
    int sz = size();
    for ( int i = 0; i < sz; i++ )
    {
        Ms.push_back(d / qs[i]);
        us.push_back(invert(Ms[i], qs[i]));
        Bs.push_back(us[i] * Ms[i]);
    }

    vint c0(sz);
    cjk = vvint(sz, c0);
    for (int i = 0; i < sz; i++)
    {
        cjk[i][i] = 0;
        for (int j = 0; j < sz; j++)
            if ( i != j )
                cjk[i][j] = invert(qs[i], qs[j]);
    }

    Y = 0;
    for ( auto x : Ms ) Y += x;
}

RnsForm Rns::split(Integer x) const
{
    assertRnsInited();
    RnsForm r;
    for (auto q : qs) r.push_back(x % q);

    if (1) if ( blend(r) != x ) never; // debug assertion
    return r;
}

Integer Rns::blend(const RnsForm & x) const
{
    int sz = size();
    auto d = dynrange();
    Integer sum(0);

    for ( int i = 0; i < sz; i++ )
    {
        auto m = x[i] * Bs[i];
        sum += m % d;
    }
    return sum % d;
}

string Rns::print() const
{
    ostringstream os;
    os << "Dynrange = " << dynrange() << '\n';
    os << "qs ="; for ( auto x : qs ) os << ' ' << x; os << '\n';
    os << "Ms ="; for ( auto x : Ms ) os << ' ' << x; os << '\n';
    os << "us ="; for ( auto x : us ) os << ' ' << x; os << '\n';
    os << "Bs ="; for ( auto x : Bs ) os << ' ' << x; os << '\n';
    os << "Y = " << Y << '\n';
    return os.str();
}


Integer Rns::EBx(RnsForm x) const
{
    Integer sum(0);
    int sz = size();
    if ( (int)x.size() != sz ) never;

    for ( int i = 0; i < sz; i++ )
    {
        sum += x[i] * Bs[i];
    }
    return sum;
}

Integer Rns::EMx(RnsForm x) const
{
    Integer sum(0);
    int sz = size();
    if ( (int)x.size() != sz ) never;

    for ( int i = 0; i < sz; i++ )
    {
        sum += (x[i] * us[i]) % qs[i] * Ms[i];
    }
    return sum;
}

vint Rns::mrs(RnsForm x) const
{
    int sz = size();
    vint v(sz), xb(sz);

    for (int k = 0; k < sz; k++)
    {
        auto q = qs[k];
        v[k] = x[k];
        for (int j = 0; j < k; j++)
        {
            auto vkj = (v[k] + q - v[j]) % q;
            v[k] = vkj * cjk[j][k] % q;
        }
        xb[k] = v[k];
    }
    return xb;
}


vint Rns::xi(RnsForm x) const
{
    Integer sum(0);
    int sz = size();
    if ( (int)x.size() != sz ) never;

    vint r(sz);

    for ( int i = 0; i < sz; i++ )
        r[i] = (x[i] * us[i]) % qs[i];

    return r;
}










string s(int x, int w)
{
    string r = std::to_string(x);
    while ( r.size() < w ) r = " " + r;
    return r;
}

string w(string r, int b = 0, int a = 0)
{
    for ( int i = 0; i < b; i++ ) r = ' ' + r;
    for ( int i = 0; i < a; i++ ) r += ' ';
    return r;
}

inline string s2(int i) { return s(i, 2); }
inline string s3(int i) { return s(i, 3); }
inline string s4(int i) { return s(i, 4); }
inline string s5(int i) { return s(i, 5); }
inline string s6(int i) { return s(i, 6); }
inline string s7(int i) { return s(i, 7); }
inline string s8(int i) { return s(i, 8); }

void cmain()
{

    if (1)
    {
        Rns rns({13, 15, 17, 19});
        //Rns rns({11, 5, 7});
        //Rns rns({3, 5, 7});
        //const int x = 11, y = 13, z = 15;
        //Rns rns({11, 13, 15});
        //const int x = 17, y = 23, z = 29;
        //Rns rns({17, 23, 29});
        //Rns rns({5, 7, 11});
        //Rns rns({7, 9});

        int sz = (int)rns.qs.size();

        auto dynrange = rns.dynrange();
        cout << "Rnse:\n" << rns.print() << '\n';

        cout << "   X| x0 x1 x2 |    XM  r   Xk  p  |  chi     |  MRS     | Headers\n";
        cout << "----------------------------------------------------------\n";
        for ( int i = 0; i < dynrange * 1 + 0; i++ )
        {
            RnsForm rf = rns.split(i);
            cout << s4(i) << w("|");
            for ( auto x : rf ) cout << s3(x);
            cout << w("|", 1);
            cout << s8(rns.EBx(rf));

            auto trueRank0 = ( rns.EBx(rf) - i );
            auto trueRank = trueRank0 / dynrange;
            if ( trueRank0 % dynrange ) never;
            cout << s3(trueRank);

            cout << s7(rns.EMx(rf));

            auto normRank0 = ( rns.EMx(rf) - i );
            auto normRank = normRank0 / dynrange;
            if ( normRank0 % dynrange ) never;
            cout << s3(normRank) << w("|", 1);

            // inverted rank
            if (0)
            {
                auto j = dynrange - i;
                if ( !i ) j = 0;
                auto normRank01 = ( rns.EMx(rns.split(j)) - j );
                auto normRank1 = normRank01 / dynrange;
                if ( normRank01 % dynrange ) never;
                cout << s3(normRank1 + normRank) << w("|", 1);
            }

            ///for ( int i = 0; i < sz; i++  ) cout << s3(rf[i]*rns.us[i] % rns.qs[i]);
            auto xis = rns.xi(rf);
            for ( auto x : xis ) cout << s3(x);

            cout << w("|", 1);
            auto mrs = rns.mrs(rf);
            for ( auto x : mrs ) cout << s3(x);
            cout << w("|", 1);

            cout << '\n';
        }

    }

    // parity maps 2D
    if (1)
    {
        cout << "\n";

        //const int x = 27, y = 67;
        //const int x = 5, y = 7;
        //const int x = 9, y = 11;
        const int x = 7, y = 9;
        Rns rns({x, y});
        cout << rns.print() << '\n';
        int sz = (int)rns.qs.size();
        int dyn = rns.dynrange();

        int t[x][y];
        for ( int i = 0; i < dyn; i++ )
        {
            RnsForm f  = rns.split(i);
            t[f[0]][f[1]] = i;
        }

        auto tableN = [&t, x, y]()
        {
            for ( int i = 0; i < x; i++ )
            {
                for ( int j = 0; j < y; j++ )
                    cout << s5(t[i][j]);
                cout << '\n';
            }
        };

        auto tableP = [&t, x, y]()
        {
            for ( int i = 0; i < x; i++ )
            {
                for ( int j = 0; j < y; j++ )
                    cout << ((t[i][j] + i + j) % 2 ? "#" : "`");
                cout << '\n';
            }
        };

        auto tableR = [&t, x, y, &rns, dyn]()
        {
            for ( int i = 0; i < x; i++ )
            {
                for ( int j = 0; j < y; j++ )
                {
                    auto a = t[i][j];
                    RnsForm rf = rns.split(a);
                    auto normRank0 = ( rns.EMx(rf) - a );
                    auto normRank = normRank0 / dyn;
                    cout << (normRank % 2 ? "#" : "`");
                }
                cout << '\n';
            }
        };

        //tableN();
        //tableP();

        for ( int i = 0; i < dyn; i++ )
        {
            vint xi  = rns.xi(rns.split(i));
            t[xi[0]][xi[1]] = i;
        }

        tableN();
        //tableP();
        tableR();
    }


    // parity maps 3D
    if (1)
    {
        cout << "\n";

        //const int x = 11, y = 31, z = 107;
        //const int x = 3, y = 5, z = 7;
        //const int x = 5, y = 7, z = 9;
        const int x = 13, y = 15, z = 17;
        Rns rns({x, y, z});
        cout << rns.print() << '\n';
        int sz = (int)rns.qs.size();
        int dyn = rns.dynrange();

        int t[x][y][z];
        for ( int i = 0; i < dyn; i++ )
        {
            RnsForm f  = rns.split(i);
            t[f[0]][f[1]][f[2]] = i;
        }

        auto tableN = [&t, x, y, z]()
        {
            for ( int i = 0; i < x; i++ )
            {
                cout << "x = " << i << '\n';
                for ( int j = 0; j < y; j++ )
                {
                    for ( int k = 0; k < z; k++ )
                        cout << s6(t[i][j][k]);
                    cout << '\n';
                }
                cout << '\n';
            }
            cout << "\n\n";
        };

        auto tableP = [&t, x, y, z]()
        {
            for ( int i = 0; i < x; i++ )
            {
                cout << "x = " << i << '\n';
                for ( int j = 0; j < y; j++ )
                {
                    for ( int k = 0; k < z; k++ )
                        cout << ((t[i][j][k] + i + j + k) % 2 ? "#" : "`");
                    cout << '\n';
                }
                cout << '\n';
            }
            cout << '\n';
        };

        //tableN();
        //tableP();

        for ( int i = 0; i < dyn; i++ )
        {
            vint xi  = rns.xi(rns.split(i));
            t[xi[0]][xi[1]][xi[2]] = i;
        }

        tableN();
        tableP();
    }


    // map 3D extended
    if (1)
    {
        cout << "\n";

        //const int x = 11, y = 31, z = 107;
        //const int x = 3, y = 5, z = 7;
        //const int x = 5, y = 7, z = 9;
        const int x = 17, y = 23, z = 29;

        const int n = 32;
        Rns rns({x, y, z});
        cout << rns.print() << '\n';
        int sz = (int)rns.qs.size();
        int dyn = rns.dynrange();

        int t[n][n][n];
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                for ( int k = 0; k < n; k++ )
                    t[i][j][k] = -1;

        for ( int i = 0; i < dyn; i++ )
        {
            RnsForm f  = rns.split(i);
            auto f0 = (f[0] * n + x / 2) / x;
            auto f1 = (f[1] * n + y / 2) / y;
            auto f2 = (f[2] * n + z / 2) / z;
            t[f0][f1][f2] = i;
        }

        auto tableN = [&t, n]()
        {
            for ( int i = 0; i < n; i++ )
            {
                cout << "x = " << i << '\n';
                for ( int j = 0; j < n; j++ )
                {
                    for ( int k = 0; k < n; k++ )
                    {
                        auto a = t[i][j][k];
                        //if( a>=0 )
                        cout << s6(a);
                    }
                    cout << '\n';
                }
                cout << '\n';
            }
            cout << "\n\n";
        };

        auto tableP = [&t, x, y, z, n, &rns, dyn]()
        {
            for ( int i = 0; i < n; i++ )
            {
                cout << "x = " << i << '\n';
                for ( int j = 0; j < n; j++ )
                {
                    for ( int k = 0; k < n; k++ )
                    {
                        auto a = t[i][j][k];
                        if ( a >= 0 )
                        {
                            vint y  = rns.xi(rns.split(a));
                            int b = (i + j + k) / (n) % 2;
                            int c = (a + y[0] + y[1] + y[2] + b);
                            int apx = ((i + j + k) % n + 1) % n;
                            if ( apx <= 2 )
                            {
                                cout << (c % 2 ? "%" : "/");
                                //cout<<"["<<(a>dyn/2?a-dyn:a)<<"]";
                            }
                            else cout << (c % 2 ? "@" : "=");
                        }
                        else cout << "`";
                    }
                    cout << '\n';
                }
                cout << '\n';
            }
            cout << '\n';
        };

        auto listV = [&t, x, y, z, n, &rns, dyn]()
        {
            for ( int i = 0; i < n; i++ )
            {
                for ( int j = 0; j < n; j++ )
                {
                    for ( int k = 0; k < n; k++ )
                    {
                        auto a = t[i][j][k];
                        if ( a >= 0 )
                        {
                            vint y  = rns.xi(rns.split(a));
                            //int b = (i+j+k)/(n)%2;
                            //int c = (a + y[0] + y[1] + y[2]+b);
                            //int apx = ((i+j+k)%n+1)%n;
                            //if( apx <=2 )
                            //{
                            //    cout << (c % 2 ? "%" : "/");
                            //    //cout<<"["<<(a>dyn/2?a-dyn:a)<<"]";
                            //}
                            //else cout << (c % 2 ? "@" : "=");
                            cout << s6(a) << w("|")
                                 << s3(y[0]) << s3(y[1]) << s3(y[2])
                                 << '\n';
                        }
                    }
                }
            }
        };

        //tableN();
        //tableP();

        for ( int i = 0; i < dyn; i++ )
        {
            vint f  = rns.xi(rns.split(i));
            auto f0 = (f[0] * n + x / 2) / x;
            auto f1 = (f[1] * n + y / 2) / y;
            auto f2 = (f[2] * n + z / 2) / z;
            t[f0][f1][f2] = i;
        }

        //tableN();
        //tableP();
        //listV();
    }
}
