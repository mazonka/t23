#pragma once

#include <vector>
#include <ostream>

using Integer = int64_t;
using vint = std::vector<Integer>;
inline int isize(vint v) { return (int)v.size(); }

template <class T>
inline std::ostream & operator<<(std::ostream & os, const std::vector<T> & v)
{
    for (auto x : v) os << ' ' << x;
    return os;
}

Integer modmul(Integer a, Integer b, Integer mod);
Integer modadd(Integer a, Integer b, Integer mod);
Integer modsub(Integer a, Integer b, Integer mod);
Integer modinv(Integer a, Integer mod);
bool coprime(Integer a, Integer b);
Integer modpow(Integer a, Integer b, Integer m);
inline double todbl(Integer a) { return double(a); }

// add counting overflows: 0 or 1
std::pair<Integer, int> modaddOver(Integer a, Integer b, Integer mod);
