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

int main()
try
{
    if (0) ntt_timing();

    if (1)
    {
        //t04_mul3_b1();
        t00_ntt();
        t01_encode();
        t02_encSk();
        t03_encPk();
        t04_mul3();
        t05_mul2_b1();
        t05_mul2();
        t06_mul1();
        t07_mul();
    }

    if (1)
    {
        t08_decomp();
        t09_hyb1();
        t10_hyb2();
    }
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



