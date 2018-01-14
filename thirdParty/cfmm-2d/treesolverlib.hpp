//#include transformation.h		//contains all the relevant transformations
//#include calcUpWardPass.cpp		//contains all the relevant transformations
//#include calcTreeDownward.cpp		//contains all the relevant transformations
//will link these files instead!

//this is a generic fmm, others can inheret from this one and redefine the two transformation functions that it uses. for now the funs will be geared towards point sources.
//~
//~ %	function [outU,outV,phi,outTree] = pointVortexVelFMM2D(vortX, vortY, vortStr, vortCoreSize, maxParticlesInBox, pValue, inIS_USING_MEX)
//~ %	------------:
//~ %	DESCRIPTION	: Algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//~ %	AUTHOR		: Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//~ %	MODIFIED	: 26.08.2008
//~ %	NOTES		: Assumes that all inputs into this function are valid and what you want.
//~ %	TODO		: Implement multi-targets/sources, dynamicP, multilevels.
//~ %	Add in input kernel function + addapt potential fields. function...add
//~ %	checks for existance of mex files and spit out warning if none are
//~ %	there...message truncate function, remove interacted list
#include "boxtree.hpp"
#include "kernellib.hpp"
#include "sys/time.h"


class FastMultipoleMethod;

typedef struct
{
    bool isUsingMultipole;
    size_t level;
    //~ size_t startBoxNum,endBoxNum;
    double minPrecision;
	BoxTree * fmmTree;
	ParticleContainer * particleData;
	PointsContainer * targetData;
	vector<IntArray*> * binarray;
	FastMultipoleMethod * thisobj;
} DownwardArguments;

typedef struct
{
    size_t level;
	BoxTree * fmmTree;
	ParticleContainer * particleData;
	vector<IntArray*> * binarray;
	FastMultipoleMethod * thisobj;
} UpwardArguments;

typedef struct
{
    size_t TOT_DIR_CNT;
    size_t DIR_CNT;
    size_t PHANT_CNT;
    size_t NEIGHB_CNT;
    bool IS_ACTIVE;
} DebugStatsStruct;

typedef struct
{
    double panelTolerance;
    bool isInfiniteSheet;
} PanelInputStruct;

void * threadedLevelDownward(void *);
void * threadedLevelUpward(void *);

class FastMultipoleMethod
{
public:
    // Con/Destructors
    FastMultipoleMethod();
    FastMultipoleMethod(bool isDebugging);

    //for binmat
    static double bincoeff(size_t n, size_t m,vector<IntArray*> &binarray);
    static void genBinCoeff(size_t pval,vector<IntArray*> &binarray);
    static double calcbincoeff(size_t n, size_t m);

    void calcTreeUpward(BoxTree &fmmTree, ParticleContainer &data,vector<IntArray*> &binarray, size_t maxThreads);
    void calcTreeDownward(BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData, vector<IntArray*> &binarray, double precision, size_t maxThreads);

    friend void * threadedLevelDownward(void *);
    friend void * threadedLevelUpward(void *);

    DebugStatsStruct dbstats;
    void * extraField;
private:

    //upward pass functions
    void normalLevelUpward(size_t ilevel, BoxTree &fmmTree, ParticleContainer &data,vector<IntArray*> &binarray);
    void levelUpward( FmmTreeCoord boxcoord, BoxTree &fmmTree, ParticleContainer &data, vector<IntArray*> &binarray);


    void calcMultipoleExpansion(FmmBoxPointer box, ParticleContainer &data);
    void translateMultipoleExpansions(FmmBoxPointer box, BoxTree &tree,vector<IntArray*> &binarray);

    //downward pass functions
    void normalLevelDownward(size_t ilevel, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &data,vector<IntArray*> &binarray, double minPrecision, bool isUsingMultipole);
    void levelDownward(FmmTreeCoord boxcoords, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData,vector<IntArray*> &binarray, double minPrecision, bool IS_USING_MULTIPOLE);

    virtual void doFinalParticleCalcs(FmmBoxPointer box, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData, double minPrecision, bool isUsingMultipole);
    ComplexDouble calcDirectPhi(ComplexDouble particlePos,double minPrecision, FmmBoxPointer box, BoxTree &fmmTree, ParticleContainer &particleData);
    ComplexDouble calcDirectFromBox(FmmBoxPointer box, ParticleContainer &IOdata, ComplexDouble evalPos, double minPrecision);

    ComplexDouble translateBltoPoint(ComplexDouble evalPos, FmmBoxPointer box);
    void shiftInfluencesToChildren(FmmBoxPointer box, BoxTree &fmmTree,vector<IntArray*> &binarray);
    void doWellSeperatedChildInteractons(FmmBoxPointer box, BoxTree &fmmTree, double minPrecision, vector<IntArray*> &binarray);
    void translateFarBLCoeff(FmmBoxPointer child, FmmBoxPointer wellSeperatedChild, double minPrecision, vector<IntArray*> &binarray);

	virtual ComplexDouble kernelVel(size_t particleIndex, ParticleContainer &IOdata, ComplexDouble evalPos, double precision)=0;
	virtual void transformFMMvel(ComplexDouble &vel)=0;
};

class LambVortexFMM: public FastMultipoleMethod
{
    public:
    static void run( double accuracy, size_t targetsPerBox, size_t maxThreads, double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets, bool withFMM);
    virtual void transformFMMvel(ComplexDouble &vel);
    virtual ComplexDouble kernelVel(size_t particleIndex, ParticleContainer &IOdata, ComplexDouble evalPos, double precision);
};

class LambVortexVortFMM: public FastMultipoleMethod
{
    public:
    static void run( double accuracy, size_t targetsPerBox, size_t maxThreads,
                     double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets);
    virtual void transformFMMvel(ComplexDouble&) {}; // This function isnt implemented, need {}
    virtual ComplexDouble kernelVel(size_t particleIndex, ParticleContainer &IOdata, ComplexDouble evalPos, double precision);
};

class VortexElementFMM: public FastMultipoleMethod
{
    public:
    static void run(double precision, size_t targetsPerBox, size_t maxThreads, double * xl, double * yl,double * strl, double * xr, double * yr,double * strr,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numsheets, size_t numtargets, double panelTolerance, double assumePointLength, bool withFMM, bool isInfiniteSheet);
    virtual void transformFMMvel(ComplexDouble &vel);
    virtual ComplexDouble kernelVel(size_t particleIndex, ParticleContainer &particleData, ComplexDouble evalPos, double precision);
};

class LambVortexMergeFMM: public FastMultipoleMethod
{
    public:
    static void run( double in_mincoreratio, double in_maxcoreratio, double in_radiusratio, double in_maxcoresizeratio, bool usingFMM, size_t targetsPerBox, size_t maxThreads, double * x, double * y,double * str, double * core, double * prev_u, double * prev_v, size_t numparticles, double * &out_x, double * &out_y,double * &out_str, double * &out_core, double * &out_prev_u, double * &out_prev_v, size_t &out_numparticles);
    virtual void transformFMMvel(ComplexDouble &) {};
    virtual ComplexDouble kernelVel(size_t , ParticleContainer &, ComplexDouble , double ) {return ComplexDouble();};

    private:
    virtual void doFinalParticleCalcs(FmmBoxPointer box, BoxTree &fmmTree, PointsContainer &, ParticleContainer &particleData, double ,  bool);
    void findMergeCandidates(FmmBoxPointer box, ParticleMergeContainer &particleData, ComplexDouble evalPos, double cuttoffRadiusSqrd, double mincore, double maxcore, IntList &mergeList);
    size_t mergeParticles(IntList &mergelist, size_t mainParticleInd, ParticleMergeContainer &particleData, double maxcoresize);
};

double timeDiff(timeval &start, timeval &end);




//inlines
inline double FastMultipoleMethod::bincoeff(size_t n, size_t m,vector<IntArray*> &binarray)
{
    //~ return binarray.at(n)->at(m);
    return binarray[n]->operator[](m);
}
