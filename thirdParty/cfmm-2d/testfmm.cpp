//==============|
//  NAME        : testfmm
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 27.07.2009
//  DESCRIPTION : Py interface to the Fmm lambVortex, vortexelement and semi-infinite vortex sheet fmm functions. a fast multipole algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//  NOTES       : .
//  TODO        :
//==============|
#include "treesolverlib.hpp"
#include <cstdlib>
using namespace std;

int main()
{
    /* initialize random seed: */
    srand ( time(NULL) );


    size_t n = 1000000;
    size_t runs = 1;
    size_t nt = n;
    size_t targetsPerBox = 15;
    size_t maxPvalue=9;
    size_t threads = 4;
    double coreSizeSqrd=0.0;

    double * x = new double[n];
    double * y = new double[n];
    double * s = new double[n];
    double * core = new double[n];

    double * ou = new double[nt];
    double * ov = new double[nt];
    double * ex = new double[nt];
    double * ey = new double[nt];

    // allocate array
    for (size_t i=0; i<n; i++)
    {
        x[i] = double(rand()) / RAND_MAX;
        y[i] = double(rand()) / RAND_MAX;
        s[i] = double(rand()) / RAND_MAX;
        core[i] = coreSizeSqrd * double(rand()) / RAND_MAX;
    }
    for (size_t i=0; i<=nt; i++)
    {
        ex[i] = x[i];//rand();
        ey[i] = y[i];//rand();
        ou[i] = 0.0;
        ov[i] = 0.0;
    }


    timeval start,stop; gettimeofday(&start, NULL);
    double timeFull;
    for (size_t i=0; i<runs; i++)
    {
        gettimeofday(&start, NULL);

        LambVortexFMM fmm;
        fmm.run( maxPvalue, targetsPerBox, threads, x, y, s, core,
                ex, ey, ou, ov, n , nt);
        gettimeofday(&stop, NULL); timeFull = timeDiff(start, stop);
        cout << i << ": full time taken " << timeFull << endl;
    }
    delete ou;
    delete ov;
    delete x;
    delete y;
    delete s;
    delete core;
    delete ex;
    delete ey;

    return 0;
}
