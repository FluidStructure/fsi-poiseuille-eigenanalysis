//==============|    
//  NAME        : lambVortexFMM.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 26.08.2008
//  DESCRIPTION : Mex interface to the Fmm lambVortex function, a fast multipole algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//  NOTES       : .
//  TODO        : 
//==============|    

#include "treesolverlib.hpp"
#include "mex.h"

// *** MATLAB interface function ***
void mexFunction( int numLeftSide, mxArray *leftSideData[], int numRightSide, const mxArray *rightSideData[] ) 
{ 
    double *u, *v;
    double *x, *y, *str, *core, *precision,*particlesPerBox,*maxTreeLevel, *maxThreads;
    size_t numparticles, numtargets;

    // Check for proper number of arguments
    if (numRightSide != 8) 
    {
       mexErrMsgTxt("Wrong number of input arguments. Is meant to be of form ( precision, maxThreads, [x], [y], [str], [core], [evalx], [evaly])"); 
    }
    if (numLeftSide != 2) 
    {
      mexErrMsgTxt("Wrong number of output arguments. Is meant to be of form [uVel, vVel]"); 
    }

    // Get matrix dimensions
    numparticles = mxGetM(rightSideData[3]);
    numtargets = mxGetM(rightSideData[7]);

    // Create MATLAB matrices for the return arguments
    leftSideData[0] = mxCreateDoubleMatrix(numtargets,1, mxREAL); 
    leftSideData[1] = mxCreateDoubleMatrix(numtargets,1, mxREAL);    

    // Assign pointers/values to the argument matrices
    u  = mxGetPr(leftSideData[0]);                  // outputs
    v  = mxGetPr(leftSideData[1]);                  // outputs

    precision = mxGetPr(rightSideData[0]);                   // inputs
    maxThreads        = mxGetPr(rightSideData[1]);
    x     = mxGetPr(rightSideData[2]);
    y    = mxGetPr(rightSideData[3]);   
    str    = mxGetPr(rightSideData[4]);
    core = mxGetPr(rightSideData[5]);
    double * ex = mxGetPr(rightSideData[6]);
    double * ey = mxGetPr(rightSideData[7]);

    // Do the actual computations using my c function
    LambVortexNaive calc;
    calc.run( *precision, *maxThreads, x,y,str,core,ex, ey, u,v,numparticles, numtargets);
    //~ double precision, size_t maxThreads,
             //~ double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets

    return;    //    Data is returned through the leftSideData pointers
}
