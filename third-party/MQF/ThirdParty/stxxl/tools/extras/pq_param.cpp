/***************************************************************************
 *  tools/extras/pq_param.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <cmath>
#include <stxxl/types>

typedef stxxl::int64 int64;

int64 D(int64 m, int64 k, int64 MaxS, int64 E, int64 B)
{
    return (m * m / 4 - (MaxS * E / (k - m)) / B);
}

using std::cout;
using std::endl;

int main()
{
    int64 IntM = 128 * 1024 * 1024;
    int64 E = 4;
    int64 B = 8 * 1024 * 1024;
    //int64 Bstep = 128 * 1024;
    int64 MaxS = (int64(128) * int64(1024 * 1024 * 1024)) / E;

    for ( ; B > 4096; B = B / 2)
    {
        int64 m = 1;
        int64 k = IntM / B;
        for ( ; m < k; ++m)
        {
            int64 c = (k - m);
            //if( D(m,k,MaxS,E,B)>= 0 && c > 10)
            if ((k - m) * m * m * B / (E * 4) >= MaxS)
            {
                cout << (k - m) * (m) * (m * B / (E * 4)) << endl;
                cout << MaxS << endl;
                cout << "D: " << D(m, k, MaxS, E, B) << endl;
                cout << "B: " << B << endl;
                cout << "c: " << c << endl;
                cout << "k: " << k << endl;
                int64 Ae = m / 2;
                int64 Ae1 = Ae + int64(sqrt((double)D(m, k, MaxS, E, B)));
                int64 Ae2 = Ae - int64(sqrt((double)D(m, k, MaxS, E, B)));
                int64 x = c * B / E;
                int64 N = x / 4096;
                cout << "Ae : " << Ae << endl;
                cout << "Ae1: " << Ae1 << endl;
                cout << "Ae2: " << Ae2 << endl;
                cout << "minAe :" << (MaxS / (x * Ae)) << endl;
                cout << "x  : " << x << endl;
                cout << "N  : " << N << endl;
                int64 Cons = x * E + B * (m / 2) + MaxS * B / (x * (m / 2));
                cout << "COns : " << Cons << endl;
                return 0;
            }
        }
    }

    cout << "No valid parameter found" << endl;
}
