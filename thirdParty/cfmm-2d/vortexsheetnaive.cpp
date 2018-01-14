//	function [uVel,vVel] = vortexelementnaive([xl], [yl], [xr], [yr], [strr], [elalX], [elalY], threads);
//	------------:
//	DESCRIPTION	: This is a generic solver for the velocity in a 2d space induced by many(or 1) semi-infinite vortex sheets.
//	AUTHOR		: Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//	MODIFIED	: 11.10.08
//	NOTES		: Positive vorticity is clockwise
//	------------:
#include "treesolverlib.hpp"
#include "mex.h"

void mexFunction( int nlhs, mxArray *leftSideData[], int nrhs, const mxArray *rightSideData[] ) { 
   unsigned int numvorts, numevals;
    
   // Check for proper number of arguments
   if (nrhs != 8) 
   {
   	mexErrMsgTxt("Wrong number of input arguments. Is meant to be of form ([xl], [yl], [xr], [yr], [strr], [elalX], [elalY], threads)"); 
   }
   if (nlhs != 2) 
   {
      mexErrMsgTxt("Wrong number of output arguments. Is meant to be of form [uVel, vVel]"); 
   }
   
   // Get matrix dimensions
   numvorts = mxGetM(rightSideData[0]);
   numevals = mxGetM(rightSideData[5]);
   
   // Create MATLAB matrices for the return arguments
   leftSideData[0] = mxCreateDoubleMatrix(numevals,1, mxREAL); 
   leftSideData[1] = mxCreateDoubleMatrix(numevals,1, mxREAL);    
    
   // Assign pointers/values to the argument matrices
   double *outUvel  = mxGetPr(leftSideData[0]);                  // outputs
   double *outVvel  = mxGetPr(leftSideData[1]);                  // outputs
   
   double *xPosl 		= mxGetPr(rightSideData[0]);                   // inputs
   double *yPosl		= mxGetPr(rightSideData[1]);
   double *xPosr 		= mxGetPr(rightSideData[2]);                   // inputs
   double *yPosr		= mxGetPr(rightSideData[3]);
   double *strengthr	= mxGetPr(rightSideData[4]);
   double *xEvalPos		= mxGetPr(rightSideData[5]);
   double *yEvalPos 	= mxGetPr(rightSideData[6]);
   double maxThreads 	= *mxGetPr(rightSideData[7]);
    
    // Do the actual computations using my c function
    VortexSheetNaive calc;
    calc.run(maxThreads, xPosl , yPosl, xPosr , yPosr, strengthr, xEvalPos, yEvalPos, outUvel, outVvel, numvorts, numevals);   
    return;	//	data has been returned through the pointers
}
