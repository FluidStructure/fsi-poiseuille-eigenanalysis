//==============|
//  NAME        : calcdownwardpass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 07.10.2008
//  DESCRIPTION : Downward pass of a treesolving merge algorithm for the ccsvm
//  NOTES       : Is not a standalone file, requires fastmultipolemethod.cpp function to call it into action.
//  TODO        : .
//==============|

#include "treesolverlib.hpp"

void LambVortexMergeFMM::doFinalParticleCalcs(FmmBoxPointer box, BoxTree &fmmTree, PointsContainer&, ParticleContainer &IOdata, double, bool)
{
    // diverge away from standard fmm algorithm now
    ParticleMergeContainer &particleData = dynamic_cast<ParticleMergeContainer &>(IOdata);
    // Cycle through all particles in this box and sum its direct influence at the point.
// Works on a single box
    //
    //
    size_t MIN_MERGE_NUMBER = 2; // itself + one other
    IntList * targetList = box->getParticleList();
    IntList mergeList = IntList();
    FmmBoxPointer nearBox;
    FmmTreeCoord nearBoxcoords;
    FmmTreeCoordList * phantomBoxList;
    size_t cndodgy, particleInd;
    cndodgy = 0;

    double mincoreratio, maxcoreratio, maxcoresizeratio, radiusratio;  
    particleData.getMergeParameters(mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio);// For each point in this box, calculate influences from multipole.

    for (IntListIterator ipart=targetList->begin(); ipart != targetList->end(); ipart++)
    {
        particleInd = *ipart;   
        // on a per target basis
        mergeList.clear();

        ComplexDouble targetPos = particleData.getPos(particleInd);
        double thiscore = particleData.getCoreSizeSqrd(particleInd);
        double mincore = mincoreratio * thiscore;
        double maxcore = maxcoreratio * thiscore;
        double maxcoresize = thiscore * maxcoresizeratio;
        double cuttoffRadius = thiscore * radiusratio;
        double cuttoffRadiusSqrd = pow(cuttoffRadius, 2);

        // ok got main target to see if merge! with others
        // For all near neighbours, including self at position 9
        for (size_t ibox=0; ibox < MAX_NEIGHBOURS; ibox++)
        {
            nearBoxcoords = box->getNeighbour(ibox);
            if (FmmBoxClass::boxExists(nearBoxcoords))
            {
                nearBox = fmmTree.getBox(nearBoxcoords);
                findMergeCandidates(nearBox, particleData, targetPos, cuttoffRadiusSqrd, mincore, maxcore, mergeList);
            }
        }
        // For all phantom boxes next
        phantomBoxList = box->getPhantomBoxList();
        for (FmmTreeCoordListIterator ibox = phantomBoxList->begin(); ibox != phantomBoxList->end(); ibox++)
        {
            nearBox = fmmTree.getBox(*ibox);
            findMergeCandidates(nearBox, particleData, targetPos, cuttoffRadiusSqrd, mincore, maxcore, mergeList);
        }


        // now we have a merge list, lets use it!! particle list modified in place...deleted vortices set to 0 strength
        if (mergeList.size() > MIN_MERGE_NUMBER)
        {
            // already done!   
            //mergeList.insert(0, particleInd);
            cndodgy += mergeParticles(mergeList, particleInd, particleData, maxcoresize);//tmpcndodgy;
        }   
        // All interactions finished for this particle...it should never get touched again!!!
    }
    
    if (cndodgy > 0)
    {
        cout << "doing " << cndodgy << " cn = maxcoresize dodgies, in treeparticlege.cpp" << endl;
    }
} // dataContainer fields modified for all particles.


void LambVortexMergeFMM::findMergeCandidates(FmmBoxPointer box, ParticleMergeContainer &particleData, ComplexDouble evalPos, double cuttoffRadiusSqrd, double mincore, double maxcore, IntList &mergeList)
{
    IntList * particleList = box->getParticleList();
    double coresize, rsqrd;
    size_t particleInd;
    // cycle through all particles and see if it meets merging criteria
    for (IntListIterator ipart=particleList->begin(); ipart != particleList->end(); ipart++)
    {
        particleInd = *ipart;   
        // only add to list if core size is within min<o<max ratios...(suggest 0.9 + 1.1)
        if (particleData.isMergeCandidate(particleInd))
        {
            coresize = particleData.getCoreSizeSqrd(particleInd);
            if ((coresize >= mincore) && (coresize <= maxcore))
            {
                // setting merge radius at some distance this current particle would split! maybe should set it on external field factors...ie merge at a radius such that it will not split next step from diffusion/travel!!.
                rsqrd = norm(particleData.getPos(particleInd) - evalPos);
                if (rsqrd < cuttoffRadiusSqrd)
                {
                    // !!!!!!USING MARKS HODGE-PODGE MERGING
                    // APPROXIMATIONS...THESE ARE possibly NOT VALID TO ROSSI ALGORITHM.
                    //if rsqrd <= R*minsigmasqrd
                    //cumulativestrength = cumulativestrength + vormat(jmerge,3);
                    //if cumulativestrength < mergingepsilon
                    mergeList.push_back(particleInd);
                }
            }
        }    // done all influences at a single calculation point
    }
}



size_t LambVortexMergeFMM::mergeParticles(IntList &mergeList, size_t mainParticleInd, ParticleMergeContainer &particleData, double maxcoresize)
{  
    // Merges to this particle are possible so doing them now
    // mergecount += 1
    // Get the strength of the new vortex, conserving 0th moment
    double stn = 0.0;
    double xn = 0.0;
    double yn = 0.0;
    double cn = 0.0;
    double rs, st;
    ComplexDouble newpos, pos;
    size_t cndodgy = 0;
    size_t particleInd;
    for (IntListIterator ipart=mergeList.begin(); ipart != mergeList.end(); ipart++)
    {
        particleInd = *ipart;   
        st = particleData.getStrength(particleInd);
        pos = particleData.getPos(particleInd);
        stn += st;
        // Get the x-position of the new vortex conserving 1st moment
        xn += st * pos.real();
        // Get the y-position of the new vortex conserving 1st moment
        yn += st * pos.imag();
    }
    yn /= stn;
    xn /= stn;

    newpos = ComplexDouble(xn, yn);

    for (IntListIterator ipart=mergeList.begin(); ipart != mergeList.end(); ipart++)
    {
        particleInd = *ipart;   
        // Get the radial distances to the old cores
        rs = norm(particleData.getPos(particleInd) - newpos); 
        // Get the core-size of the new vortex
        // cn = sqrt(sum(strg.*(4*(cors.^2) + Rs))/(4*stn));         // This is from Rossi
        cn += particleData.getStrength(particleInd) * (pow(particleData.getCoreSizeSqrd(particleInd),2) + rs);          // This is modified for my core (no 4)  !!!!!!! ABSOLUTE VALUE!!!!!!
    }
    
    cn = sqrt(cn/stn);
    // Impose a limit on the maximum core size that can be created
    if (cn > maxcoresize)
    {                                                       
        //!!!!!!!! CAREFUL HERE TOO !!!!!!!!!!!!
        cn = maxcoresize;
        cndodgy += 1;
    }

    // modifying vortex in current position
    particleData.setMergeCandidate(mainParticleInd, false);
    particleData.setPos(mainParticleInd, newpos);
    particleData.setStrength(mainParticleInd, stn);
    particleData.setCoresize(mainParticleInd, cn);

    for (IntListIterator ipart=mergeList.begin(); ipart != mergeList.end(); ipart++)
    {
        particleInd = *ipart;   
        // deleting merged vortices except for the original which should be at the back
        if (particleInd != mainParticleInd)
        {
            particleData.setMergeCandidate(particleInd, false);
            particleData.setStrength(particleInd, 0.0);
        }
    }

    return cndodgy;
}
