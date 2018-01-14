//	function [uVel,vVel] = vortexelementnaive([xl], [yl], [strl],[xr], [yr], [strr], [elalX], [elalY], precision, threads);
//	------------:
//	DESCRIPTION	: This is a generic solver for the velocity in a 2d space induced by many(or 1) finite vortex sheets.
//	AUTHOR		: Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//	MODIFIED	: 29.09.08
//	NOTES		: Positive vorticity is clockwise
//	TODO		: .
//	------------:
#include "treesolverlib.hpp"
#include "mex.h"

void mexFunction( int nlhs, mxArray *leftSideData[], int nrhs, const mxArray *rightSideData[] ) { 
   unsigned int numvorts, numevals;
    
   // Check for proper number of arguments
   if (nrhs != 10) 
   {
   	mexErrMsgTxt("Wrong number of input arguments. Is meant to be of form ([xl], [yl], [strl],[xr], [yr], [strr], [elalX], [elalY], precision, threads)"); 
   }
   if (nlhs != 2) 
   {
      mexErrMsgTxt("Wrong number of output arguments. Is meant to be of form [uVel, vVel]"); 
   }
   
   // Get matrix dimensions
   numvorts = mxGetM(rightSideData[0]);
   numevals = mxGetM(rightSideData[6]);
   
   // Create MATLAB matrices for the return arguments
   leftSideData[0] = mxCreateDoubleMatrix(numevals,1, mxREAL); 
   leftSideData[1] = mxCreateDoubleMatrix(numevals,1, mxREAL);    
    
   // Assign pointers/values to the argument matrices
   double *outUvel  = mxGetPr(leftSideData[0]);                  // outputs
   double *outVvel  = mxGetPr(leftSideData[1]);                  // outputs
   
   double *xPosl 		= mxGetPr(rightSideData[0]);                   // inputs
   double *yPosl		= mxGetPr(rightSideData[1]);
   double *strengthl	= mxGetPr(rightSideData[2]);
    double *xPosr 		= mxGetPr(rightSideData[3]);                   // inputs
   double *yPosr		= mxGetPr(rightSideData[4]);
   double *strengthr	= mxGetPr(rightSideData[5]);
   double *xEvalPos		= mxGetPr(rightSideData[6]);
   double *yEvalPos 	= mxGetPr(rightSideData[7]);
   double precision 	= *mxGetPr(rightSideData[8]);
   double maxThreads 	= *mxGetPr(rightSideData[9]);
    // Do the actual computations using my c function
    VortexElementNaive calc;
    calc.run( precision, maxThreads, xPosl , yPosl, strengthl, xPosr , yPosr, strengthr, xEvalPos, yEvalPos, outUvel, outVvel, numvorts, numevals);   

    return;	//	data has been returned through the pointers
}
