//==============|
//  NAME        : calcdownwardpass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 07.10.2008
//  DESCRIPTION : Downward pass of a fast multipole algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//  NOTES       : Is not a standalone file, requires fastmultipolemethod.cpp function to call it into action.
//  TODO        : .
//==============|

//  transformationfunction is function that will transfer the phi,u,v into the one for the correct application.  outU = vVel/(2*pi());outV = -uVel/(2*pi());    //   this turns now returns in the correct format for a clockwise vorticity field
//     corefunction: is the function to apply to the "direct" calculations, ie, the naive ones because they are too close to use the fmm for.

//    ------------:
//    NAME:         calcTreeDownward
//    PURPOS:     Performs the entire downward pass of the FMM
//    IMPORTS:     .
//                .
//                .
//    PRE-CONDS:     outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//    POST-CONDS: outPhi, outU, outV contain valid data
//    ------------:
#include "treesolverlib.hpp"
pthread_mutex_t downwardLevelMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t dbstatsMutex = PTHREAD_MUTEX_INITIALIZER;
size_t downwardBoxCounter = 0; //used for the threads checking out their own box...must be reset every level

void FastMultipoleMethod::calcTreeDownward(BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData, vector<IntArray*> &binarray, double precision, size_t maxThreads)
{
    size_t THREADS=maxThreads;
    size_t minboxthreshold = 1;

    FmmBoxPointer box;
    FmmTreeCoord boxcoords;

    boxcoords.level = 0;
    boxcoords.index = 1;

    box = fmmTree.getBox(boxcoords);

    bool IS_USING_MULTIPOLE = false;

    if (box->getPvalue() > 0)
    {
        IS_USING_MULTIPOLE = true;
    }

    double minPrecision = precision;

    for (size_t ilevel = 0; ilevel <= fmmTree.getMaxLevel(); ilevel++)
    {
        // Start at level 1 because multipole method doesnt work at level 1 (points no longer are well seperated)
        // otherwise the children wont have a bl set
        downwardBoxCounter = 0;     // ooooooops that fucked me!
        size_t totalBoxes = fmmTree.getNumBoxes(ilevel);

        if ((THREADS > 1) && (totalBoxes > minboxthreshold))
	    {
            pthread_t threadID[THREADS];
            // sorting out input arguments
            DownwardArguments threadArgs;
            threadArgs.isUsingMultipole = IS_USING_MULTIPOLE;
            threadArgs.level = ilevel;
            threadArgs.minPrecision = minPrecision;
            threadArgs.fmmTree = &fmmTree;
            threadArgs.particleData = &particleData;
            threadArgs.targetData = &targetData;
            threadArgs.binarray = &binarray;
            threadArgs.thisobj = this;

            //~ DownwardArguments threadArgs[THREADS];

            for(size_t i=0; i < THREADS; i++)
            {
                //~ threadArgs[i].startBoxNum = 1 + (int)( i * totalBoxes/THREADS);
                //~ threadArgs[i].endBoxNum = (int)((i+1) * totalBoxes/THREADS);
                //~ threadArgs[i].level = ilevel;
                //~ threadArgs[i].minPrecision = minPrecision;
                //~ threadArgs[i].fmmTree = &fmmTree;
                //~ threadArgs[i].particleData = &particleData;
                //~ threadArgs[i].binarray = &binarray;

                //~ pthread_create( &threadID[i], NULL, threadedLevelDownward, (void *) &threadArgs[i] );
                pthread_create( &threadID[i], NULL, threadedLevelDownward, (void *) &threadArgs );
            }

            //make sure threads have joined before moving on
            for(size_t i=0; i < THREADS; i++)
            {
                pthread_join( threadID[i], NULL);
            }
        }
        else
        {
            normalLevelDownward(ilevel,fmmTree, targetData, particleData,binarray,minPrecision, IS_USING_MULTIPOLE);
        }

    } // No more levels, downward pass completed outPhi, outU, outV contain their valid particleData.

    return;
}

void * threadedLevelDownward(void * inarguments)
{
    FmmTreeCoord boxcoords;
    DownwardArguments * arguments = static_cast<DownwardArguments *>(inarguments);

    bool IS_USING_MULTIPOLE = arguments->isUsingMultipole;
    size_t ilevel = arguments->level;

    double minPrecision = arguments->minPrecision;
    ParticleContainer &particleData = *arguments->particleData;
    PointsContainer &targetData = *arguments->targetData;
	BoxTree & fmmTree= *arguments->fmmTree;
	vector<IntArray*> binarray = *arguments->binarray;
	FastMultipoleMethod * thisobj = arguments->thisobj;

    boxcoords.level = ilevel;

	size_t totalBoxes = fmmTree.getNumBoxes(ilevel);
    pthread_mutex_lock( &downwardLevelMutex);
    while (downwardBoxCounter < totalBoxes)
    {
        // cycle through all boxes on this level
        // checkout this box so no other threads try and do the same one!
        downwardBoxCounter++;
        boxcoords.index = downwardBoxCounter;
        pthread_mutex_unlock( &downwardLevelMutex);
        //other threads are now able to continue the loop

        // Doing the exciting bit now
        thisobj->levelDownward(boxcoords, fmmTree, targetData, particleData, binarray, minPrecision, IS_USING_MULTIPOLE);

        pthread_mutex_lock( &downwardLevelMutex);
    } // No more boxes on this level
    pthread_mutex_unlock( &downwardLevelMutex);

    return NULL;
}

void FastMultipoleMethod::normalLevelDownward(size_t ilevel, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData,vector<IntArray*> &binarray, double minPrecision, bool IS_USING_MULTIPOLE)
{
    FmmTreeCoord boxcoords;

    boxcoords.level = ilevel;
    for (size_t ibox=1; ibox <= fmmTree.getNumBoxes(ilevel); ibox++)
    {
        // cycle through all boxes on this level
        boxcoords.index = ibox;
        levelDownward(boxcoords, fmmTree, targetData, particleData, binarray, minPrecision, IS_USING_MULTIPOLE);
    } // No more boxes on this level
    return;
}

void FastMultipoleMethod::levelDownward(FmmTreeCoord boxcoords, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData,vector<IntArray*> &binarray, double minPrecision, bool IS_USING_MULTIPOLE)
{
    FmmBoxPointer box = fmmTree.getBox(boxcoords);
    if  (box->holdsTargets())       // Dont care about downward motion if theres no targets
    {
        if (box->isLowestBox())
        {
            // No more children, the particles are kept in this box.
            doFinalParticleCalcs(box, fmmTree, targetData, particleData, minPrecision, IS_USING_MULTIPOLE);    // Calculate final phi,u,v at the targets (Greengard STEP 5/6/7/8)
        }
        else if (IS_USING_MULTIPOLE)
        {
            // only do each childs translation if that child itself has targets
            // Box is a parent that has children, must filter the multipole at its centre to the children and keep recursing.
            shiftInfluencesToChildren(box, fmmTree, binarray);                                    // Shift translate blCoeff from parents to children. (Greengard STEP 4)

            // Regardless if this box is lowest level and contains vortices, or it contains children boxes, always do the well seperated (long distance) multipole traslations to its center first. This is the crux of the FMM
            doWellSeperatedChildInteractons(box, fmmTree, minPrecision, binarray);    // shifting far neighbour aKCoeffs from the upward pass onto itself. (Greengard STEP 3)
        }
    } // No more actions left to do inside this box.
}








//============================================================================
//=============== NEAR-FIELD/DIRECT CALCULATIONS =============================
//============================================================================

//    ------------:
//    NAME:         doFinalParticleCalcs
//    PURPOS:     Performs the entire downward pass of the FMM
//    IMPORTS:     .
//                .
//                .
//    PRE-CONDS:     outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//    POST-CONDS: outPhi, outU, outV contain valid particleData
//    ------------:

void FastMultipoleMethod::doFinalParticleCalcs(FmmBoxPointer box, BoxTree &fmmTree, PointsContainer &targetData, ParticleContainer &particleData, double minPrecision,  bool IS_USING_MULTIPOLE)
{
    // Works on a single box

            // First shift the bl onto the targets,
        //    In the boxes blCoeff it captures all of the influences
        //    from the well seperated points at that level

            // Influence from boxes multipole back onto the targets
            // [phi(vortInd),uVel(vortInd),vVel(vortInd)] = calcWellSeparatedPhi(boxBl,zVel,(evalX + 1i*evalY),pValue);

            //shiftWellSeparatedPhi(box, targetData);
            //shifting influences from the multipole at centre of this box!

    ComplexDouble targetPos,vel;
    IntList * targetList;

    targetList = box->getTargetList();

    // For each point in this box, calculate influences from multipole.
    for (IntListIterator ipart=targetList->begin(); ipart != targetList->end(); ipart++)
    {

        // on a per target basis
        targetPos = targetData.getPos(*ipart);
        if (IS_USING_MULTIPOLE)
        {
            vel = translateBltoPoint(targetPos, box); // returns vel as if it was a simple source
           // Need to apply our transform from logZ into whatever we want
            transformFMMvel(vel);
        }
        else
        {
            vel=0;
        }

            //    The remaining influences that have not been accounted for
            //    are obtained through the standard direct Phi method...this
            //    is what the FMM is about, reducing the number of these calcs!
            //    have kept this in same format to take Dr Mark Pitman's vectorized direct vort code.

            // Find influence from particles in current box first
        vel += calcDirectPhi(targetPos, minPrecision, box, fmmTree, particleData);

        targetData.setVel(*ipart,vel);
        // All interactions finished for this particle...it should never get touched again!!!
    } // dataContainer fields modified for all particles.

}




//  ------------:
//  NAME        : translateBltoPoint
//  PURPOSE     : Performs the entire downward pass of the FMM
//  IMPORTS     : .
//
//  PRE-CONDS   : outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//  POST-CONDS  : outPhi, outU, outV contain valid data
//  ------------:
ComplexDouble FastMultipoleMethod::translateBltoPoint(ComplexDouble evalPos, FmmBoxPointer box)
{
    ComplexDouble zRel, bl, velocity;
    size_t pValue;

    pValue = box->getPvalue();
    zRel = evalPos -  box->getPos();
    //    taken from lemma 2.??? to shift the local multipole back onto a local point
    // velocity is a simple dphi/dz
    velocity=0;
    for(size_t l = 1; l <= pValue; l++)
    {
        bl = box->getBlCoeff(l);
        velocity += (double)l * bl * pow(zRel,(l-1));
    }

    //return the conjugate of the veloicty calculated
    return conj(velocity);
}

//  ------------:
//  NAME        : calcDirectPhi
//  PURPOSE     : Shifts blcoeffs and phantom boxes from parent to children
//  IMPORTS     : .
//                .
//                .
//  PRE-CONDS   : outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//  POST-CONDS  : outPhi, outU, outV contain valid data
//  ------------:
ComplexDouble FastMultipoleMethod::calcDirectPhi(ComplexDouble targetPos, double minPrecision, FmmBoxPointer box, BoxTree & fmmTree, ParticleContainer &particleData)
{
    bool debuggingIsActive = this->dbstats.IS_ACTIVE;
    ComplexDouble velocity(0,0);
    FmmBoxPointer evalbox;
    FmmTreeCoord evalboxcoords;
    FmmTreeCoordList * phantomBoxList;
    // For all near neighbours, including self at position 9
    for (size_t ibox=0; ibox < MAX_NEIGHBOURS; ibox++)
    {
        evalboxcoords = box->getNeighbour(ibox);
        if (FmmBoxClass::boxExists(evalboxcoords))
        {
            evalbox = fmmTree.getBox(evalboxcoords);
            velocity += calcDirectFromBox(evalbox, particleData, targetPos, minPrecision);

            if (debuggingIsActive)
            {
                pthread_mutex_lock( &dbstatsMutex);
                this->dbstats.DIR_CNT += evalbox->getNumParticles();
                pthread_mutex_unlock( &dbstatsMutex);
            }
        }
    }

    // For all phantom boxes next
    phantomBoxList = box->getPhantomBoxList();
    for (FmmTreeCoordListIterator ibox = phantomBoxList->begin(); ibox != phantomBoxList->end(); ibox++)
    {
        evalbox = fmmTree.getBox(*ibox);
        velocity += calcDirectFromBox(evalbox, particleData, targetPos, minPrecision);
        if (debuggingIsActive)
        {
            pthread_mutex_lock( &dbstatsMutex);
            this->dbstats.PHANT_CNT += evalbox->getNumParticles();
            pthread_mutex_unlock( &dbstatsMutex);
        }
    }

    return velocity;
}


//  ------------:
//  NAME        : calcDirectFromBox
//  PURPOS      : Takes in a whole box, and a location its direct influence is calculated upon, needs dataContainer to know what particles are in its box and where.
//  IMPORTS     : box: Pointer to a particle containing box.
//              : particleData: real data for particles in the box.
//              : evalPos: where the direct interaction is to be calculated upon.
//  PRE-CONDS   : box has a correct list of particleindexes matching to the particle container.
//  POST-CONDS  : the velocity will be the direct influences of all particles in this box at the eval point.
//  ------------:
ComplexDouble FastMultipoleMethod::calcDirectFromBox(FmmBoxPointer box, ParticleContainer &particleData, ComplexDouble evalPos, double minPrecision)
{
    // Cycle through all particles in this box and sum its direct influence at the point.
    //double particleStr, coreSize, precision;
    //ComplexDouble vel , particlePos;
    IntList * particleList;
    //bool debuggingIsActive = this->dbstats.IS_ACTIVE;

    particleList = box->getParticleList();
    ComplexDouble vel = ComplexDouble(0,0);
    for (IntListIterator ipart=particleList->begin(); ipart != particleList->end(); ipart++)
    {
        // on a per particle basis
        /*
        if debuggingIsActive
        {
            pthread_mutex_lock( &dbstatsMutex);
            this->dbstats.TOT_DIR_CNT += 1;
            pthread_mutex_unlock( &dbstatsMutex);
        }
        */
        vel += kernelVel(*ipart, particleData, evalPos, minPrecision);
    }    // done all influences at a single calculation point

   return vel;
}








//============================================================================
//=============== PARENT->CHILD TRANSLATIONS =================================
//============================================================================

//  ------------:
//  NAME        : shiftInfluencesToChildren
//  PURPOS      : Shifts blcoeffs and phantom boxes from parent to children
//  IMPORTS     : .
//  PRE-CONDS   :     outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//  POST-CONDS  : outPhi, outU, outV contain valid data
//  TODO        : Bug in lemma 2.5 missing a z^l
//  ------------:
void FastMultipoleMethod::shiftInfluencesToChildren(FmmBoxPointer box, BoxTree &fmmTree, vector<IntArray*> &binarray)
{
    ComplexDouble boxPos, childPos, childBl, parentBk, zRel;
    FmmTreeCoord childCoord;
    FmmTreeCoordList * phantomBoxList;
    FmmBoxPointer childBox;
    size_t pValue;

    boxPos = box->getPos();
    pValue = box->getPvalue();
    phantomBoxList = box->getPhantomBoxList();

    for(size_t i = 0; i < MAX_CHILDREN; i++)
    {
        // Cycle through all 4 children
        childCoord = box->getChild(i);
        if (FmmBoxClass::boxExists(childCoord))
        {
            childBox = fmmTree.getBox(childCoord); // !!!!!previous seg fault for not detecting this one!!!
            if (childBox->holdsTargets())   // dont care children with no targets
            {
                // This child actually has targets in it an hence got created/exists.
                if ( box->getNumPhantomBoxes() > 0 )
                {
                    childBox->setPhantomBoxlist(phantomBoxList);
                } // Passed phantom particle information down a generation

                // translateLocalBlCoeff(childPos, boxPos, boxBlCoeff, pValue, binCoeffMat3, childBlCoeff)
                // Using Greengard lemma 2.5
                childPos = childBox->getPos();
                zRel = boxPos-childPos;
                // Shifts blcoeffs from parent to child, Using Greengard lemma 2.5
                for (size_t l = 0; l <= pValue; l++)
                {
                    childBl = ComplexDouble(0);
                    for (size_t k = l; k <= pValue; k++)
                    {
                        parentBk = box->getBlCoeff(k);
                        childBl += parentBk * bincoeff(k,l, binarray) * pow(-zRel,(k-l)); // !!!!!!is this missing a * zl?
                    }

                    childBox->setBlCoeff(l,childBl);  // blCoeff now set for this child box. first time too so save can neglect previous bl values
                }
            }
        }
    } // No more children to account for, box interactions over
    return;
}






//============================================================================
//=============== WELL SEPERATED MULTIPOLE TRANSLATIONS ======================
//============================================================================

//    ------------:
//    NAME:         doWellSeperatedChildInteractons
//    PURPOS:     Performs the entire downward pass of the FMM
//    IMPORTS:     .
//                .
//                .
//    PRE-CONDS:     outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//    POST-CONDS: outPhi, outU, outV contain valid data
//    ------------:
void FastMultipoleMethod::doWellSeperatedChildInteractons(FmmBoxPointer box, BoxTree &fmmTree, double minPrecision, vector<IntArray*> &binarray)
{
    // This couldnt have been done the step before because these boxes were not in the interaction list. Any neighbours that are too close this time will be dealt with at a lower level(higher resolution.
    // step 3 is the long distance translation of multipoles in the interaction
    // cant have more than this rogue particles
    //this field is to be passed down from generation to generation until reach the bottom!this is
    bool debuggingIsActive = this->dbstats.IS_ACTIVE;

    ComplexDouble childPos,neighbourChildPos;

    FmmBoxPointer neighbourbox, child, neighbourChild;
    FmmTreeCoord neighbourcoords, childcoord, neighbchildcoords;


    for (size_t ineighb = 0; ineighb < MAX_NEIGHBOURS; ineighb++)
    {
        neighbourcoords = box->getNeighbour(ineighb);
        if (FmmBoxClass::boxExists(neighbourcoords))
        {
            neighbourbox = fmmTree.getBox(neighbourcoords);
            if (neighbourbox->holdsParticles())
            {
                // Neighbour exists and has particles, so it will either be a phantom box, or it will have farfield children that will have particles...but some may not.
                if (neighbourbox->isLowestBox())
                {
                    // Neighbour has no children to do a far field interaction with this boxes children, need to add this box to the lonely uncle list, so can remember to do a direct interaction for the particles in this childrens box later.
                    for (size_t ichild = 0; ichild < MAX_CHILDREN; ichild++)
                    {
                        // Cycle through all 4 children
                        childcoord = box->getChild(ichild);
                        if (FmmBoxClass::boxExists(childcoord))
                        {
                            // Child exists, but this parent neighbour has no children hence its a phantom box
                            child = fmmTree.getBox(childcoord);
                            if (child->holdsTargets())
                            {
                                // may be able to save some memory by not doing this when theres no targets, hence will never be used.
                                child->addPhantomBox(neighbourcoords);
                            }
                        }
                    }
                }
                else
                {
                    // Is a children bearing neighbour so must cycle through its children and do the actual interactions where necessary
                    for (size_t ineighbchild = 0; ineighbchild < MAX_CHILDREN; ineighbchild++)
                    {
                        neighbchildcoords = neighbourbox->getChild(ineighbchild);
                        if (FmmBoxClass::boxExists(neighbchildcoords))
                        {
                            neighbourChild = fmmTree.getBox(neighbchildcoords);
                            if (neighbourbox->holdsParticles())
                            {
                                // no point doing far field if there is no particles, it means there no AK either.
                                neighbourChildPos = neighbourChild->getPos();


                                // doing it for all the children of the main box in question!
                                // DONT NEED TO DO THIS IF THE CHILD HAS NO TARGETS WITHIN IT!!
                                for (size_t ichild = 0; ichild < MAX_CHILDREN; ichild++)
                                {
                                    // Cycle through all 4 children
                                    childcoord = box->getChild(ichild);
                                    if (FmmBoxClass::boxExists(childcoord))
                                    {
                                        child = fmmTree.getBox(childcoord);
                                        if (child->holdsTargets())
                                        {
                                            if (child->isWellSeperated(neighbourChild))
                                            {
                                                // Cousin exists and is well seperated.
                                                // translating well seperated aks into bls
                                                // important to remember dont need a BL coefficient if there is no targets in the box, as it will never get used!!!                         this
                                                if (debuggingIsActive)
                                                {
                                                    pthread_mutex_lock( &dbstatsMutex);
                                                    this->dbstats.NEIGHB_CNT += 1;
                                                    pthread_mutex_unlock( &dbstatsMutex);
                                                }
                                                translateFarBLCoeff(child, neighbourChild, minPrecision, binarray);          // Translates + increments the BlCoeff from the far field.
                                            }
                                        }
                                    }
                                } // All of this neighbours-children accounted for far-field wise
                            }
                        }
                    } // All of the neighbours accounted for this current childs Bl coefficient is valid.
                }
            }
        }
    } // No more children to account for, box interactions over
}

//  ------------:
//  NAME:         translateFarBLCoeff
//  PURPOS:     Performs the entire downward pass of the FMM
//  IMPORTS:     .
//              .
//              .
//  PRE-CONDS:     outPhi, outU, outV contain pointers to empty arrays of valid length. fmmTree has valid akCoeffs after an upward pass
//  POST-CONDS: outPhi, outU, outV contain valid data
//  ------------:
void FastMultipoleMethod::translateFarBLCoeff(FmmBoxPointer child, FmmBoxPointer wellSeperatedChild, double minPrecision, vector<IntArray*> &binarray)
{
    ComplexDouble zRel, bl, childPos, wellSeperatedChildPos, farAzero, farAk, kseries;
    double cVal, childBoxRsize;
    size_t dynPvalue;

    childPos = child->getPos();
    childBoxRsize = child->getRsize();
    wellSeperatedChildPos = wellSeperatedChild->getPos();

    // Using dynamicP optimizations
    cVal = abs(wellSeperatedChildPos-childPos)/childBoxRsize - 1.0;
    dynPvalue = static_cast<size_t>(1 + static_cast<int>((-log(1.00001*minPrecision))/log(cVal)));

    zRel = wellSeperatedChildPos-childPos;

    //    taken from lemma 1.??? to shift the well seperated multipoles onto a
    //    local point, therefore converting Ak to a bl (from taylor series to a polynomial series)

    //    Calcing Blzero first
    farAzero = wellSeperatedChild->getAkCoeff(0);
    bl = farAzero * log(-1.0 * zRel);

    for(size_t k = 1; k <= dynPvalue; k++)
    {
        farAk = wellSeperatedChild->getAkCoeff(k);
        bl += ( farAk / pow(zRel,k) ) * pow(-1.0, (int)k);
        // Maybe bad here...have a divide by zRel
    }
    bl += child->getBlCoeff(0); // preserving original state!

    child->setBlCoeff(0, bl);

    //    Now calculating remainingbl for l>0
    for(size_t l = 1; l <= dynPvalue; l++)
    {
        bl = -1.0 * farAzero / ((double)l * pow(zRel,l) );    //    some mixed mode arith going on here int*double

        kseries = ComplexDouble (0,0);
        for (size_t k = 1; k <= dynPvalue; k++)
        {
            farAk = wellSeperatedChild->getAkCoeff(k);
            kseries += farAk * bincoeff(l+k-1,k-1, binarray) * pow(-1.0, (int)k) / pow(zRel,k) ;
        }    //    Taylor series now summed up for this current l

        bl += kseries / pow(zRel,l);
        bl += child->getBlCoeff(l); // preserving original state!

        child->setBlCoeff(l, bl);
    }
   return;    //    data has been set
}
