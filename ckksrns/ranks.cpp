

#include "ccrun.h"

#include "integer.cpp"

Integer mod;

Integer power(Integer x, int pw)
{
    if ( x == 0 ) return 0;
    if ( pw == 0 ) return 1;

    if ( pw < 0)
    {
        x = modinv(x, mod);
        pw = -pw;
    }

    Integer r {1};
    while (pw--) r = modmul(r, x, mod);
    return r;
}

Integer y(Integer x, int pw)
{
    auto a = power(x, pw);
    //return a;
    return x * a / mod;
}

void cmain()
{
    mod = 257;

    for ( Integer i = 0; i < mod; i++ )
    {
        auto j = i?mod - i:0;
        cout << i
             << '\t' << y(i, -1)
             //<< '\t' << (mod%(i?i:1))
             << '\t' << y(i, 1)
             << '\t' << y(i, 2)
             << '\t' << y(i, 3)
             << '\t' << y(i, 4)
             << '\t' << y(i, 5)
             << '\t' << y(i, 6)
             << '\t' << y(i, 7)
             << '\n';
    }

}
