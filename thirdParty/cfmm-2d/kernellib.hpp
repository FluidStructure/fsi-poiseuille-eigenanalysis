//==============|    
//  NAME        : kernellibclass.hpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 30.08.2008
//  DESCRIPTION : Holds the kernel functions for many different particle interactions.
//  NOTES       : Suitable for use with the fmm or as a standalone function.
//  TODO        : .
//==============|   

#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

typedef complex<double> ComplexDouble;

class KernelLib
{
    public:
    static ComplexDouble lambVortex(ComplexDouble particlePos, double particleStr, double coreSizeSqrd, ComplexDouble evalPos, double precision);
    static ComplexDouble lambVortexVort(ComplexDouble particlePos, double particleStr, double coreSizeSqrd, ComplexDouble evalPos, double precision);
    static ComplexDouble vortexElement(ComplexDouble posL, double strL, ComplexDouble posR, double strR, ComplexDouble evalPos, double panel_tolerance, bool isInfiniteSheet);
};

inline double sign(double a)
{
    return static_cast<double>((a >= 0.0) - (a < 0.0));
}
