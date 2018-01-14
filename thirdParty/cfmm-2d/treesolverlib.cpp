//	function [outU,outV,phi,outTree] = pointVortexVelFMM2D(particleX, particleY, particleStr, vortCoreSize, maxParticlesInBox, pValue, inIS_USING_MEX)
//	------------:
//	DESCRIPTION	: Algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//	AUTHOR		: Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//	MODIFIED	: 16.06.2009
//	NOTES		: Assumes that all inputs into this function are valid and what you want.

#include "treesolverlib.hpp"

using namespace std;
bool DEBUG_ON = false;
#ifdef _DEBUG_MODE
    DEBUG_ON = true;
#endif

//===========================================================================
//=============== generic fast multipole method =============================
//===========================================================================
FastMultipoleMethod::FastMultipoleMethod()
{
    this->dbstats.IS_ACTIVE = false;
    this->dbstats.TOT_DIR_CNT = 0;
    this->dbstats.DIR_CNT = 0;
    this->dbstats.PHANT_CNT = 0;
    this->dbstats.NEIGHB_CNT = 0;
    this->extraField = NULL;
}

FastMultipoleMethod::FastMultipoleMethod(bool isDebugging)
{
    this->dbstats.IS_ACTIVE = isDebugging;
    this->dbstats.TOT_DIR_CNT = 0;
    this->dbstats.DIR_CNT = 0;
    this->dbstats.PHANT_CNT = 0;
    this->dbstats.NEIGHB_CNT = 0;
    this->extraField = NULL;
}

void FastMultipoleMethod::genBinCoeff(size_t pval,vector<IntArray*> &binarray)
{
    binarray.clear();
    binarray.resize(2*pval);
    IntArray* column;
    for (size_t n=0; n < (2*pval); n++)
    {
        // add this new vect to vector
        column = new IntArray(pval+1,0);
        binarray.at(n) = column;
        for (size_t m=0; m < (pval+1); m++)
        {
            column->at(m) = calcbincoeff(n,m);
        }
    }
    return;
}

//FROM http://www.brpreiss.com/books/opus4/html/page467.html
double FastMultipoleMethod::calcbincoeff(size_t n, size_t m)
{
   size_t b[n+1];
   b[0]=1;
   for(size_t i=1; i<=n; ++i)
   {
      b[i] = 1;
        for(size_t j=(i-1); j>0 ; --j)
        {
            b[j] += b[j-1];
        }
   }
   return b[m];
}









//===========================================================================
//=============== LambVortexFMM =================================
//===========================================================================
void LambVortexFMM::run( double accuracy, size_t targetsPerBox, size_t maxThreads, double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets, bool withFMM)
{
    timeval tstart,tstop;
    gettimeofday(&tstart, NULL);
    timeval start,stop;

    // place particles in a valid container
    ParticleContainer particleData = ParticleContainer(x,y,str,core,numparticles);
    PointsContainer targetData = PointsContainer(evalx, evaly,evalu,evalv,numtargets);

    // vars needed to be set before going into fmm:
    size_t maxTreeLevel, minTreeLevel, pvalue;
    double precision; 
    vector<IntArray*> binarray;
    if (!withFMM)
    {
        pvalue = 0; // no farfield effects needed;
        precision = 0.0;
        targetsPerBox = 0; // Keep splitting to whatever level I say;

        maxTreeLevel = 0;    // all in 1 big box...as naive as it gets!
        // DONT FORGET WANT A MAXTREELEVEL with alreast 1box/thread!
        if (maxThreads > 1)
        {
            maxTreeLevel = 1;    // cant be higher than level 1 otherwise some of the far field interactions will be missed!
        }
        minTreeLevel = maxTreeLevel; // build a tree to this exact level.
    }
    else
    {

    //determine precision and pvalues from the inputted accuracy field...can be in either a pval or actual precision
    if (accuracy >= 1)
    {
        // input was a valid pValue was
        pvalue = static_cast<size_t> (accuracy);
        precision = pow((4.0/(sqrt(2.0))-1.0),-accuracy);
    }
	else
	{
	    // input was a decimal precision value
		precision = accuracy;
		pvalue = static_cast<size_t> (1 + static_cast<int> (-log(1.00001 * precision) / log((4.0/(sqrt(2.0))-1.0))));
    }

    // maxbox size is dictated by lamb core size, sqrt because input is coresqrd
    double rCuttoff = sqrt(particleData.getMaxSize()) * sqrt(-log(precision));
    if ( rCuttoff > 0 )
    {
        maxTreeLevel = BoxTree::boxSizeToTreeLevel(rCuttoff, particleData, targetData);
    }
    else
    {
        maxTreeLevel = 50; // arbitrary maximum
    }

    minTreeLevel = 2;    // no well seperated for less than this so fmm would be waste of time

    genBinCoeff(pvalue,binarray);
    }

    //http://cboard.cprogramming.com/showthread.php?t=101085
    gettimeofday(&start, NULL);
    BoxTree tree = BoxTree();
    tree.buildAdaptive(particleData, targetData, pvalue, targetsPerBox, minTreeLevel, maxTreeLevel, maxThreads);
    gettimeofday(&stop, NULL);
    double timeBuild = timeDiff(start, stop);

    LambVortexFMM fmm;
    fmm.dbstats.IS_ACTIVE = DEBUG_ON;
    gettimeofday(&start, NULL);
    if (withFMM)
    {
        fmm.calcTreeUpward(tree,particleData,binarray, maxThreads);
    }
    gettimeofday(&stop, NULL);
    double timeUp = timeDiff(start, stop);

    gettimeofday(&start, NULL);
    fmm.calcTreeDownward(tree, targetData, particleData, binarray, precision, maxThreads);
    gettimeofday(&stop, NULL);
    double timeDown = timeDiff(start, stop);

    if (DEBUG_ON)
    {
        cout << "\t...in treesolverlib.cpp (~LN167)";
        cout << "Time taken (build-up-down) (s) , ";
        cout << timeBuild << ", " << timeUp << ", " << timeDown << endl;
        cout << "max treeleveel " << tree.getMaxLevel() << endl;
        cout << "Tree level sizes is: (maxlevel:" << tree.getMaxLevel() << ") ";
        cout << tree.levelSizesToString() << endl;
        cout << "\t...out treesolverlib.cpp (~LN167)" << endl;
        if (false)
        {
            cout << "operation count for: total direct, phantom direct, neighb";
            cout << "-neighb";
            size_t cntlist[] = {fmm.dbstats.TOT_DIR_CNT, fmm.dbstats.DIR_CNT , fmm.dbstats.PHANT_CNT, fmm.dbstats.NEIGHB_CNT};
            for(size_t i=0; i < 4; i++)
            {
                cout << ", " << cntlist[i];
            }
        }
        cout << endl;
        gettimeofday(&tstop, NULL);
        cout << "All function time " << timeDiff(tstart, tstop) << endl;
    }

    //deleting binarray
    for(size_t i=0; i < binarray.size(); i++)
    {
        delete binarray.at(i);
    }
    return;
}

//  ------------:
//  NAME        : transformFMMvel
//  PURPOSE     : is function that will transfer the phi,u,v into the one for the correct application.
//                 currently correct format for a clockwise vorticity field.
//  IMPORTS     : .
//  PRE-CONDS   : vel is the velocity from the multipole of a source/sink
//  POST-CONDS  : vel is converted to the equivalent vortex model
//  TODO        : a validation/test function tests multipole to the naive by building a 2particle far - tree with p=high
//  ------------:
inline void LambVortexFMM::transformFMMvel(ComplexDouble &vel)
{
    // !!!!!removed /2pi here!!!! 4 weeks of frustration!!!
    double u = vel.imag();
    double v = -vel.real();
    vel = ComplexDouble(u,v);
}

//  ------------:
//  NAME        : kernelVel
//  PURPOS      : Takes in a whole box, and a location its direct influence is calculated upon, needs dataContainer to know what particles are in its box and where.
//  IMPORTS     : box: Pointer to a particle containing box.
//              : IOdata: real data for particles in the box.
//              : evalPos: where the direct interaction is to be calculated upon.
//  PRE-CONDS   : box has a correct list of particleindexes matching to the particle container.
//  POST-CONDS  : the velocity will be the direct influences of all particles in this box at the eval point.
//  TODO        : precisions!
//  ------------:
inline ComplexDouble LambVortexFMM::kernelVel(size_t particleIndex, ParticleContainer &particleData, ComplexDouble evalPos, double precision)
{
    ComplexDouble pos = particleData.getPos(particleIndex);
    double str = particleData.getStrength(particleIndex);
    double core = particleData.getCoreSizeSqrd(particleIndex);
	// Cycle through all particles in this box and sum its direct influence at the point.
    return KernelLib::lambVortex(pos, str , core, evalPos, precision);
}




//===========================================================================
//=============== LambVortexVortFMM =================================
//===========================================================================
void LambVortexVortFMM::run( double accuracy, size_t targetsPerBox, size_t maxThreads, double * x, double * y,double * str, double * core,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numparticles, size_t numtargets)
{
	// place particles in a valid container
    ParticleContainer particleData = ParticleContainer(x,y,str,core,numparticles);
    PointsContainer targetData = PointsContainer(evalx, evaly,evalu, evalv ,numtargets);

    //determine precision and pvalues from the inputted accuracy field...can be in either a pval or actual precision

    double precision;
    size_t pvalue;
    pvalue = 0; // Dont care about upward pass or any far field effects, only near field is significant. Still need precision however
    if (accuracy >= 1)
    {
        // input was mistakenly an error or pvalue
        precision = pow((4.0/(sqrt(2.0))-1.0),-accuracy);
    }
	else
	{
	    // input was a decimal precision value
		precision = accuracy;
    }

    // maxbox size is dictated by lamb core size ....does this return max size or max size squared....i think it will be squared??? so sqrt
    double rCuttoff = sqrt(particleData.getMaxSize()) * sqrt(-log(precision));

    size_t maxTreeLevel;
    if ( rCuttoff > 0 )
    {
        maxTreeLevel = BoxTree::boxSizeToTreeLevel(rCuttoff, particleData, targetData);
    }
    else
    {
        maxTreeLevel = 50; // arbitrary maximum
    }

    size_t minTreeLevel = 2;    // no well seperated for less than this so fmm would be waste of time

    /////////////////for debugging!!
    // NOTE, THE MULTILEVELLING HAS BEEN DISABLED FOR THIS CODE AS IT HAD BUGS WHEN THERE MINTREE != MAXTREE LEVEL
    minTreeLevel = maxTreeLevel; // build a tree to this exact level.
   // THIS IS A NAIVE-ish CODE
    ///////////////////

    vector<IntArray*> binarray;
    //cout << "Time taken (build-up-down) (s) , ";

    timeval start,stop;
    //~ //http://cboard.cprogramming.com/showthread.php?t=101085
    gettimeofday(&start, NULL);
    BoxTree tree = BoxTree();
    tree.buildAdaptive(particleData, targetData, pvalue, targetsPerBox, minTreeLevel, maxTreeLevel, maxThreads);
    gettimeofday(&stop, NULL);
    // cout << timeDiff(start, stop) << ", ";
    LambVortexVortFMM fmm;
    fmm.dbstats.IS_ACTIVE = DEBUG_ON;

    // dont need upward pass for this, as farfield vorticity is approximated as 0

    gettimeofday(&start, NULL);
    fmm.calcTreeDownward(tree, targetData, particleData, binarray, precision, maxThreads);
    gettimeofday(&stop, NULL);
    //cout << timeDiff(start, stop) << ", ";

    return;
}

//  ------------:
//  NAME        : kernelVel
//  PURPOS      : Takes in a whole box, and a location its direct influence is calculated upon, needs dataContainer to know what particles are in its box and where.
//  IMPORTS     : box: Pointer to a particle containing box.
//              : IOdata: real data for particles in the box.
//              : evalPos: where the direct interaction is to be calculated upon.
//  PRE-CONDS   : box has a correct list of particleindexes matching to the particle container.
//  POST-CONDS  : the velocity will be the direct influences of all particles in this box at the eval point.
//  TODO        : precisions!
//  ------------:
inline ComplexDouble LambVortexVortFMM::kernelVel(size_t particleIndex, ParticleContainer &particleData, ComplexDouble evalPos, double precision)
{
    ComplexDouble pos = particleData.getPos(particleIndex);
    double str = particleData.getStrength(particleIndex);
    double core = particleData.getCoreSizeSqrd(particleIndex);
	// Cycle through all particles in this box and sum its direct influence at the point.
    return KernelLib::lambVortexVort(pos, str , core, evalPos, precision);
}







//===========================================================================
//=============== VortexElementFMM =================================
//===========================================================================
void VortexElementFMM::run( double accuracy, size_t targetsPerBox, size_t maxThreads, double * xl, double * yl,double * strl, double * xr, double * yr,double * strr,  double * evalx, double * evaly,  double * evalu, double * evalv, size_t numsheets, size_t numtargets, double panelTolerance, double assumePointLength, bool withFMM, bool isInfiniteSheet)
{ 
    // place elements and points in a valid container
    SheetContainer sheetData = SheetContainer(xl,yl,strl, xr,yr, strr, numsheets);
    PointsContainer targetData = PointsContainer(evalx,evaly,evalu,evalv,numtargets);
    PanelInputStruct extraInput;
    extraInput.panelTolerance = panelTolerance;
    extraInput.isInfiniteSheet = isInfiniteSheet;

    // vars needed to be set before going into fmm:
    size_t maxTreeLevel, minTreeLevel, pvalue;
    double precision; 
    vector<IntArray*> binarray;

    if (!withFMM || isInfiniteSheet)
    {
        pvalue = 0; // no farfield effects needed;
        precision = 0.0;
        targetsPerBox = 0; // Keep splitting to whatever level I say;

        maxTreeLevel = 0;    // all in 1 big box...as naive as it gets!
        // DONT FORGET WANT A MAXTREELEVEL with atleast 1box/thread!
        if (maxThreads > 1)
        {
            maxTreeLevel = 1;    // cant be higher than level 1 otherwise some of the far field interactions will be missed!
        }
        minTreeLevel = maxTreeLevel; // build a tree to this exact level.
    }
    else
    {
        // determine precision and pvalues from the inputted accuracy field...can be in either a pval or actual precision

        if (accuracy >= 1)
        {
            // input was a valid pValue was
            pvalue = static_cast<size_t> (accuracy);
            precision = pow((4.0/(sqrt(2.0))-1.0),-accuracy);
        }
        else
        {
            // input was a decimal precision value
            precision = accuracy;
            pvalue = static_cast<size_t> (1 + static_cast<int> (-log(1.00001 * precision) / log((4.0/(sqrt(2.0))-1.0))));
        }

        // maxbox size is dictated by the maximum sheet size

        // A FUDGE IN USE BELOW!
        double rCuttoff = assumePointLength * sheetData.getMaxLength();     // FUDGE TO SAY THAT A PANEL APPEARS AS A POINT SOURCE WHEN 4*L AWAY...THIS IS CORRECT TO 0.4%...USE 9L TO MAKE IT 0.1%
        // A FUDGE IN USE ABOVE!

        if ( rCuttoff > 0 )
        {
            maxTreeLevel = BoxTree::boxSizeToTreeLevel(rCuttoff, sheetData, targetData);
        }
        else
        {
            maxTreeLevel = 30; // arbitrary maximum
        }

        minTreeLevel = 2;    // no well seperated for less than this so fmm would be waste of time

        genBinCoeff(pvalue,binarray);
    }

    // build the quad tree
    BoxTree tree = BoxTree();
    tree.buildAdaptive(sheetData, targetData, pvalue, targetsPerBox, minTreeLevel, maxTreeLevel, maxThreads);
    if (DEBUG_ON)
    {
        cout << "max treeleveel " << tree.getMaxLevel() << endl;
    }
    VortexElementFMM fmm;
    fmm.extraField = static_cast<void *>(&extraInput);
    fmm.dbstats.IS_ACTIVE = DEBUG_ON;
    if (withFMM)
    {
        // run upward pass, treating sheets as points
        fmm.calcTreeUpward(tree, sheetData, binarray, maxThreads);
    }
    // run the calcs on the tree
    fmm.calcTreeDownward(tree, targetData, sheetData, binarray, precision, maxThreads);

    if (withFMM)
    {
        //deleting binarray
        for(size_t i=0; i < binarray.size(); i++)
        {
            delete binarray.at(i);
        }
    }
    return;
}

//  ------------:
//  NAME        : transformFMMvel
//  PURPOSE     : is function that will transfer the phi,u,v into the one for the correct application.
//                 currently correct format for a clockwise vorticity field.
//  IMPORTS     : .
//  PRE-CONDS   : vel is the velocity from the multipole of a source/sink
//  POST-CONDS  : vel is converted to the equivalent vortex model
//  TODO        : a validation/test function tests multipole to the naive by building a 2particle far - tree with p=high
//  ------------:
inline void VortexElementFMM::transformFMMvel(ComplexDouble &vel)
{
    // !!!!!removed /2pi here!!!! 4 weeks of frustration!!!
    double u = vel.imag();
    double v = -vel.real();
    vel = ComplexDouble(u,v);
}

//  ------------:
//  NAME        : kernelVel
//  PURPOS      : Takes in a whole box, and a location its direct influence is calculated upon, needs dataContainer to know what particles are in its box and where.
//  IMPORTS     : box: Pointer to a particle containing box.
//              : IOdata: real data for particles in the box.
//              : evalPos: where the direct interaction is to be calculated upon.
//  PRE-CONDS   : box has a correct list of particleindexes matching to the particle container.
//  POST-CONDS  : the velocity will be the direct influences of all particles in this box at the eval point.
//  TODO        : precisions!
//  ------------:
inline ComplexDouble VortexElementFMM::kernelVel(size_t particleIndex, ParticleContainer &IOdata, ComplexDouble evalPos, double)
{
    // get extra info out of a god pointer first
    PanelInputStruct * extraField = static_cast<PanelInputStruct *>(this->extraField);
    double panelTolerance = extraField->panelTolerance;
    bool isInfiniteSheet = extraField->isInfiniteSheet;
    SheetContainer &sheetData = dynamic_cast<SheetContainer &>(IOdata);
    ComplexDouble posl = sheetData.getPosL(particleIndex);
    ComplexDouble posr = sheetData.getPosR(particleIndex);
    double strl = sheetData.getStrengthL(particleIndex);
    double strr = sheetData.getStrengthR(particleIndex);
    // Cycle through all particles in this box and sum its direct influence at the point.
    return KernelLib::vortexElement(posl, strl, posr, strr, evalPos, panelTolerance, isInfiniteSheet);
}



//===========================================================================
//=============== LambVortexMergeFMM =================================
//===========================================================================
    void LambVortexMergeFMM::run( double in_mincoreratio, double in_maxcoreratio, double in_radiusratio, double in_maxcoresizeratio, bool usingFMM, size_t targetsPerBox, size_t maxThreads, double * x, double * y,double * str, double * core, double * prev_u, double * prev_v, size_t numparticles, double * &out_x, double * &out_y,double * &out_str, double * &out_core, double * &out_prev_u, double * &out_prev_v, size_t &out_numparticles)

{
    // place particles in a valid container
    ParticleMergeContainer particleData = ParticleMergeContainer(x,y,str,core,numparticles, in_mincoreratio, in_maxcoreratio, in_radiusratio, in_maxcoresizeratio);
    //PointsContainer targetData = PointsContainer((NULL,NULL,NULL,NULL,0); // bad, causes tree to blow out!
    PointsContainer targetData = PointsContainer(x,y,NULL,NULL,numparticles);

    double precision = 0;
    size_t pvalue = 0; // Dont care about upward pass or any far field effects, only near field is significant. Still need precision however

    // maxbox size is dictated by maximum lamb core size ....does this return max size or max size squared....i think it will be squared???
    // represents the biggest particles having a particle within its adjacent box
    double rCuttoff = particleData.getMaxSize() * in_radiusratio;

    size_t maxTreeLevel;
    if ( rCuttoff > 0 )
    {
        maxTreeLevel = BoxTree::boxSizeToTreeLevel(rCuttoff, particleData, targetData);
    }
    else
    {
        maxTreeLevel = 50; // arbitrary maximum
    }

    if (!(usingFMM) && (maxTreeLevel > 1))
    {
        maxTreeLevel = 1;
    }

    size_t minTreeLevel = maxTreeLevel; // build a tree to this exact level.

    vector<IntArray*> binarray;

    //int b = cin.get();

    BoxTree tree = BoxTree();
    tree.buildAdaptive(particleData, targetData, pvalue, targetsPerBox, minTreeLevel, maxTreeLevel, maxThreads);

    cout << "max treeleveel " << tree.getMaxLevel() << endl;

    //int a = cin.get();
    LambVortexMergeFMM fmm;
    // not thread safe yet!
    maxThreads = 1;
    fmm.calcTreeDownward(tree, targetData, particleData, binarray, precision, maxThreads);

    // now done merging...time to delete out the old cruft!! (str = 0.0)
    size_t * tmpsort = new size_t[numparticles];
    // could easily write an inplace code for more memory/speed efficiency
    size_t numkeepers = 0;

    for (size_t i=0; i<numparticles; i++)
    {
        if (particleData.getStrength(i) != 0.0)
        {
            tmpsort[numkeepers] = i;
            numkeepers++;
        }    
    }

    // creating new out space
    out_x = new double[numkeepers];
    out_y = new double[numkeepers];
    out_str = new double[numkeepers]; 
    out_core = new double[numkeepers];
    size_t realInd;
    bool withPrev = ((prev_u != NULL) && (prev_v != NULL));
    if (withPrev)
    {   
        out_prev_u = new double[numkeepers];
        out_prev_v = new double[numkeepers];
    }
    for (size_t i=0; i<numkeepers; i++)
    {
        realInd = tmpsort[i];

        out_x[i] = particleData.getPos(realInd).real();
        out_y[i] = particleData.getPos(realInd).imag();
        out_str[i] = particleData.getStrength(realInd);
        out_core[i] = particleData.getCoreSizeSqrd(realInd);
        if (withPrev)
        {
            if (particleData.isMergeCandidate(realInd))
            {
                // not part of a merge, keep existing prev_u/v
                out_prev_u[i] = prev_u[realInd];
                out_prev_v[i] = prev_v[realInd];
            }
            else
            {   
                // Means it was a part of a merge!
                out_prev_u[i] = 0.0;
                out_prev_v[i] = 0.0;
            }
        }
            
    }

    out_numparticles = numkeepers;

    delete tmpsort;
    
    return;
}

double timeDiff(timeval &start, timeval &end)
{
    double ms= (end.tv_usec - start.tv_usec);
    double s = (end.tv_sec - start.tv_sec);
    return (s+(ms/1000000));
}
