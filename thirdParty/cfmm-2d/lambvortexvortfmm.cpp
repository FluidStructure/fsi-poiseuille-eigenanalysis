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
    double *x, *y, *ex, *ey, *str, *core, *pval,*particlesPerBox,*maxTreeLevel, *maxThreads;
    size_t numparticles, numtargets;
    mxArray * tmpv;

    // Check for proper number of arguments
    if (numRightSide != 9) 
    {
       mexErrMsgTxt("Wrong number of input arguments. Is meant to be of form ( pvalue, particlesPerBox,  maxThreads, [x], [y], [str], [core], [evalx], [evaly])"); 
    }
    if (numLeftSide != 1) 
    {
      mexErrMsgTxt("Wrong number of output arguments. Is meant to be of form [scalarVorticity]"); 
    }

    // Get matrix dimensions
    numparticles = mxGetM(rightSideData[3]);
    numtargets = mxGetM(rightSideData[8]);

    // Create MATLAB matrices for the return arguments
    leftSideData[0] = mxCreateDoubleMatrix(numtargets,1, mxREAL); 
    tmpv            = mxCreateDoubleMatrix(numtargets,1, mxREAL);    

    // Assign pointers/values to the argument matrices
    u  = mxGetPr(leftSideData[0]);                  // outputs
    v  = mxGetPr(tmpv);                  // outputs

    pval = mxGetPr(rightSideData[0]);                   // inputs
    particlesPerBox        = mxGetPr(rightSideData[1]);
    maxThreads     = mxGetPr(rightSideData[2]);
    x     = mxGetPr(rightSideData[3]);
    y    = mxGetPr(rightSideData[4]);   
    str    = mxGetPr(rightSideData[5]);
    core = mxGetPr(rightSideData[6]);
    ex    = mxGetPr(rightSideData[7]);
    ey = mxGetPr(rightSideData[8]);

    // Do the actual computations using my c function
    LambVortexVortFMM fmm;
    //static void run( double accuracy, size_t particlesPerBox, size_t maxThreads, double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets);
    fmm.run( *pval, *particlesPerBox, *maxThreads, x,y,str,core, ex, ey, u,v,numparticles, numtargets);

    return;    //    Data is returned through the leftSideData pointers
}
