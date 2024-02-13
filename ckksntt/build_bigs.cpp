#include "ccrun.h"
#include "ccrut.h"

void run1(string p)
{
    sys("rm -rf main.exe o");
    sys("mv o." + p + " o");
    sys("mv main." + p + ".exe main.exe");
    sys("m -j 4 MPIR=" + p);
    sys("mv o o." + p);
    sys("mv main.exe main." + p + ".exe");
}

void cmain()
{
    run1("i");
    run1("1");
}
