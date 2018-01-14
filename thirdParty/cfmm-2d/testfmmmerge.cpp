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


    size_t numparticles = 20000;
    size_t targetsPerBox = 15;
    size_t maxThreads = 4;
    double coreSizeSqrd=0.1;
    double in_mincoreratio = 0.9;
    double in_maxcoreratio = 1.1;
    double in_radiusratio = 1.0;
    double in_maxcoresizeratio = 2.0;
    bool usingFMM = true;
    double * x = new double[numparticles];
    double * y = new double[numparticles];
    double * str = new double[numparticles];
    double * core = new double[numparticles];
    double * pu = NULL; //new double[numparticles];
    double * pv = NULL; //new double[numparticles];
    
    double * out_x, * out_y, * out_str, * out_core, *out_pu, *out_pv; 
    size_t out_numparticles;

    // allocate array
    for (size_t i=0; i<numparticles; i++)
    {
        x[i] = double(rand()) / RAND_MAX;
        y[i] = double(rand()) / RAND_MAX;
        str[i] = double(rand()) / RAND_MAX;
        core[i] = coreSizeSqrd * double(rand()) / RAND_MAX;
    }

    timeval start,stop; gettimeofday(&start, NULL);

    LambVortexMergeFMM fmm;
    fmm.run(in_mincoreratio, in_maxcoreratio, in_radiusratio, in_maxcoresizeratio, usingFMM, targetsPerBox, maxThreads, x, y,str, core, pu,pv, numparticles, out_x, out_y,out_str, out_core, out_pu, out_pv, out_numparticles);
    gettimeofday(&stop, NULL); double timeFull = timeDiff(start, stop);
    cout << "full time taken " << timeFull << endl;

    cout << "numparts before " << numparticles << " after " << out_numparticles << endl;
    cout << "trying to access memory..." << out_x[0] << endl << out_x[out_numparticles-1] << endl;


    delete x;
    delete y;
    delete str;
    delete core;
    delete out_x;
    delete out_y;
    delete out_str;
    delete out_core;

    return 0;
}
