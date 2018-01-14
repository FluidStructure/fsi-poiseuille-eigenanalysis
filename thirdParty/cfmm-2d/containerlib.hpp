//==============|    
//  NAME        : containerlib.hpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 01.09.2008
//  DESCRIPTION : Headers to ParticleContainer.cpp, a class for the state of a point/particle in 2D-space. 
//  NOTES       : 
//  TODO        : 
//==============|    
//~ #ifndef SSTREAM_H
//~ #define SSTREAM_H
#include <sstream>
//~ #endif
//~ #ifndef IOSTREAM_H
//~ #define IOSTREAM_H
#include <iostream>
//~ #endif
//~ #ifndef COMPLEX_H
//~ #define COMPLEX_H
#include <complex>
//~ #endif
using namespace std;
typedef complex<double> ComplexDouble;

class PointsContainer 
{
public:
    PointsContainer(double * inX, double *inY, double * inU, double *inV, size_t numParticles);
    PointsContainer(){};
    virtual ComplexDouble getPos(size_t partnum) { return ComplexDouble(x[partnum],y[partnum]); };   // use the last defined getPos function
    ComplexDouble getVel(size_t partnum) { return ComplexDouble(u[partnum],v[partnum]); };
    size_t getNumParticles() { return numParticles; };
    double getDomainSize();
    ComplexDouble getDomainCentre();
    void getExtremities(double &minX, double &maxX, double &minY, double &maxY );

    void setNumParticles(size_t in_numParticles) { numParticles = in_numParticles; };
    void setPos(size_t partnum, ComplexDouble inPos);
    void setVel(size_t partnum, ComplexDouble inVel);
    const string toString(size_t partnum);
    const string toString();
    
    static void getDomain(PointsContainer &particleData, PointsContainer &targetData, double &size, ComplexDouble &centre);
    void calcExtremities(double &minX, double &maxX, double &minY, double &maxY, size_t numParticles);
    void setExtremities(double minXa, double maxXa, double minYa, double maxYa);

private:
    double *x, *y, *u, *v;
    double minX, maxX, minY, maxY;
    size_t numParticles;
};

class ParticleContainer: public PointsContainer
{
public:
    ParticleContainer(double * inX, double *inY, double * inStrength, double * inCoreSizeSqrd, size_t numParticles);
    ParticleContainer(){};
    double getCoreSizeSqrd(size_t partnum) { return coreSizeSqrd[partnum]; };
    void setCoresize(size_t partnum, double core) { coreSizeSqrd[partnum] = core; };

    virtual double getStrength(size_t partnum) { return strength[partnum]; };
    virtual void setStrength(size_t partnum, double str) { strength[partnum] = str; };
    
    double getMaxSize();

    const string toString(size_t partnum);
    const string toString();
private:
    double *strength, *coreSizeSqrd;
};


class ParticleMergeContainer: public ParticleContainer
{
public:
    ParticleMergeContainer(double * inX, double *inY, double * inStrength, double * inCoreSizeSqrd, size_t inNumParticles, double inMincoreratio, double inMaxcoreratio, double inRadiusratio, double inMaxcoresizeratio);
    ParticleMergeContainer(){};
    ~ParticleMergeContainer();
    
    void getMergeParameters(double &mincoreratio, double &maxcoreratio, double &radiusratio, double &maxcoresizeratio);// For each point in this box, calculate influences from multipole.
    bool isMergeCandidate(size_t partnum) { return mergeCandidate[partnum]; };
    void setMergeCandidate(size_t partnum, bool decision) { mergeCandidate[partnum] = decision; };

private:
    bool * mergeCandidate;
    double mincoreratio, maxcoreratio, radiusratio , maxcoresizeratio;
};

class SheetContainer: public ParticleContainer
{
public:
    SheetContainer(double * inXl, double *inYl, double * inStrengthl, double * inXr, double *inYr, double * inStrengthr, size_t numSheets);
    double getStrengthL(size_t partnum) { return leftSide.getStrength(partnum); };
    double getStrengthR(size_t partnum) { return rightSide.getStrength(partnum); };
    ComplexDouble getPosL(size_t partnum) { return leftSide.getPos(partnum); };
    ComplexDouble getPosR(size_t partnum) { return rightSide.getPos(partnum); };
    double getLength(size_t partnum) { return abs(getPosR(partnum) - getPosL(partnum)); };
    
    virtual double getStrength(size_t partnum) { return getLength(partnum) * (getStrengthL(partnum) + getStrengthR(partnum)) / (2.0); };
    virtual ComplexDouble getPos(size_t partnum) { return ((getPosL(partnum) + getPosR(partnum)) / 2.0); };   // use the last defined getPos function
   
    double getMaxLength();    

    const string toString(size_t partnum);
    const string toString();
private:
    ParticleContainer leftSide, rightSide;
};

extern ostream & operator<<(ostream &os, PointsContainer &data);
extern ostream & operator<<(ostream &os, ParticleContainer &data);
extern ostream & operator<<(ostream &os, SheetContainer &data);

