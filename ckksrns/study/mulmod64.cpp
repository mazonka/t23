#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#define Assert assert
#define uint64 uint64_t

/*
 * Calculate (x * y) % m, where x and y in [0, 2^64), m in [1, 2^64).
 *
 * If x or y is greater than 2^32, improved interleaved modular
 * multiplication algorithm is used to avoid overflow.
 */
static uint64 modular_multiplicate(uint64 x, uint64 y, const uint64 m)
{
    int     i, bits;
    uint64      r = 0;

    Assert(1 <= m);

    /* Because of (x * y) % m = (x % m * y % m) % m */
    if (x >= m)
        x %= m;
    if (y >= m)
        y %= m;

    /* Return the trivial result. */
    if (x == 0 || y == 0 || m == 1)
        return 0;

    /* Return the result if (x * y) can be multiplicated without overflow. */
    if ((x | y) < (0xffffffff))
        return (x * y) % m;

    /* To reduce the for loop in the algorithm below. */
    if (x < y)
    {
        uint64 tmp = x;
        x = y;
        y = tmp;
    }

    /* Interleaved modular multiplication algorithm
     *
     *   D.N. Amanor, et al, "Efficient hardware architecture for
     *    modular multiplication on FPGAs", in Field Programmable
     *    Logic and Apllications, 2005. International Conference on,
     *    Aug 2005, pp. 539-542.
     *
     * This algorithm is usually used in the field of digital circuit
     * design.
     *
     * Input: X, Y, M; 0 <= X, Y <= M;
     * Output: R = X *  Y mod M;
     * bits: number of bits of Y
     * Y[i]: i th bit of Y
     *
     * 1. R = 0;
     * 2. for (i = bits - 1; i >= 0; i--) {
     * 3.   R = 2 * R;
     * 4.   if (Y[i] == 0x1)
     * 5.       R += X;
     * 6.   if (R >= M) R -= M;
     * 7.   if (R >= M) R -= M;
     *   }
     *
     * In Steps 3 and 5, overflow should be avoided.
     * Steps 6 and 7 can be instead of a modular operation (R %= M).
     */

    bits = 64;

    for (i = bits - 1; i >= 0; i--)
    {
        if (r > 0x7fffffffffffffff)
            /* To avoid overflow, transform from (2 * r) to
             * (2 * r) % m, and further transform to
             * mathematically equivalent form shown below:
             */
            r = m - ((m - r) << 1);
        else
            r <<= 1;

        if ((y >> i) & 0x1)
        {
            if (r > 0xffffffffffffffffull - x)
                /* To calculate (r + x) without overflow, transform to (r + x) % m,
                         * and transform to mathematically equivalent form (r + x - m).
                 */
                r += x - m;
            else
                r += x;
        }

        r %= m;
    }

    return r;
}

int main(int argc, char ** argv)
{

    if (argc != 4)
    {
        printf("Syntax Error:\n\tUsage:%s A B N\n", argv[0]);
        return -1;
    }

    uint64_t a  = strtoull(argv[1], NULL, 10);
    uint64_t b  = strtoull(argv[2], NULL, 10);
    uint64_t n  = strtoull(argv[3], NULL, 10);

    uint64_t r =  modular_multiplicate(a, b, n);
    printf("(%llu * %llu) %% %llu = %llu\n", a, b, n, modular_multiplicate(a, b, n));
    return 0;
}
