//==============|
//  NAME        : boxtreebuild.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 20.10.2008
//  DESCRIPTION : seperate implementaiton file for the FmmTree building process.
//==============|

#include "boxtree.hpp"
#include <pthread.h>

size_t numparticles;

pthread_mutex_t splitmutex = PTHREAD_MUTEX_INITIALIZER;
FmmTreeCoordListIterator splitBoxIterator;

pthread_mutex_t addBoxMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t futureSplitMutex = PTHREAD_MUTEX_INITIALIZER;

void BoxTree::buildAdaptive(PointsContainer &particleData, PointsContainer &targetData, size_t inpValue, size_t targetsPerBox, size_t minTreeLevel, size_t maxTreeLevel,size_t maxThreads)
{
    bool IS_USING_MULTIPOLE = false;
    if (inpValue > 0)
    {
        IS_USING_MULTIPOLE = true;
    }

    size_t THREADS=maxThreads;
    bool IS_DYNAMIC_R = IS_USING_MULTIPOLE;

    size_t level, initializesize;
	FmmTreeCoordList splitBoxes, futureSplitBoxes;
	FmmTreeCoord firstbox;
	numparticles = particleData.getNumParticles();
	firstbox = initializeFirstBox(particleData, targetData, inpValue);
    // sets up all needed to seed the recursive process
    splitBoxes.push_back(firstbox);
    // now add this box to the tree to say it needs to be split!
    level = 0;
    //	Giving the recursive algorithm a starting point at level 0
    //	now have vorts seperated into an initial quadrant can start dynamically allocating vortexIndex into boxes
    //	want to go through each box and find which child it sits in relative to this box. if child index = 0 then grab a new index from admin mat and then increment it

    while ((level < maxTreeLevel) && (!(splitBoxes.empty())))
	{
	    // Boxes are awaiting splitting
	    initializesize = MAX_CHILDREN * splitBoxes.size();
        // best guess for amount of memory needed, cant be any more than 4 children per box to split!
	    this->initializeLevel(level+1, initializesize);
        // avoids poor performance

		// atleast one box has too many particles in it, and we arent at maximum level yet
		if (THREADS > 1 && (level>0))
	    {
            pthread_t threadID[THREADS];
            SplitBoxArguments threadArgs;
            splitBoxIterator = splitBoxes.begin();
            // this resets the global iterator
            threadArgs.splitBoxes = &splitBoxes;
            threadArgs.particleData = &particleData;
            threadArgs.targetData = &targetData;
            threadArgs.IS_DYNAMIC_R = IS_DYNAMIC_R;
            threadArgs.thisobj = this;

            for(size_t i=0; i < THREADS; i++)
            {
                pthread_create( &threadID[i], NULL, threadedBoxSplit, (void *) &threadArgs );
            }

            // make sure threads have joined before moving on
            for(size_t i=0; i < THREADS; i++)
            {
                pthread_join( threadID[i], NULL);
            }
        }
        else
        {
            normalBoxSplit(splitBoxes, particleData, targetData, IS_DYNAMIC_R);
        }

		// Need to do one more loop through to get neighbours. DO NOT PUT THIS IN THE LOOP ABOVE!
		// Also Cycles through all children in a box and adds it to future split boxes if necessary

        if (THREADS > 1 && (level>0))
	    {
            pthread_t threadID[THREADS];
            PostSplitArguments threadArgs;
            splitBoxIterator = splitBoxes.begin();
            // this resets the global iterator
            threadArgs.splitBoxes = &splitBoxes;
            threadArgs.futureSplitBoxes = &futureSplitBoxes;

            threadArgs.targetsPerBox = targetsPerBox;
            threadArgs.MIN_TREE_LEVEL = minTreeLevel;
            threadArgs.thisobj = this;

            for(size_t i=0; i < THREADS; i++)
            {
                pthread_create( &threadID[i], NULL, threadedPostSplit, (void *) &threadArgs );
            }
            // make sure threads have joined before moving on
            for(size_t i=0; i < THREADS; i++)
            {
                pthread_join( threadID[i], NULL);
            }
        }
        else
        {
            normalPostSplit(splitBoxes, futureSplitBoxes, targetsPerBox, minTreeLevel);
        }

        splitBoxes.clear();
		splitBoxes = futureSplitBoxes;
		futureSplitBoxes.clear();
		level++;
	} // tree finished!

    bool debug_on = false;
    if (debug_on)
    {
        if (level >= maxTreeLevel)
        {
            cout << "in boxtreebuild: " << endl;
            cout << "clipping at level " << level << endl;
            cout << "in boxtreebuild: " << endl;
        }
    }
}
//  ------------:
//  NAME        : .
//  PURPOSE     : Sets up neighbours for all children.
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:


void * threadedBoxSplit(void * inarguments)
{
    FmmTreeCoord currentBoxCoords;

    SplitBoxArguments * arguments = (SplitBoxArguments *) inarguments;

	FmmTreeCoordList &splitBoxes = *arguments->splitBoxes;
	PointsContainer &particleData = *arguments->particleData;
    PointsContainer &targetData = *arguments->targetData;
   	bool IS_DYNAMIC_R = arguments->IS_DYNAMIC_R;
	BoxTree * thisobj = arguments->thisobj;

    // prevent possible race condition to the while loop
    pthread_mutex_lock( &splitmutex);
    while (splitBoxIterator != splitBoxes.end())
    {
        // lets check us out a box in a thread friendly manner
        currentBoxCoords = *splitBoxIterator;
        splitBoxIterator++;
        pthread_mutex_unlock( &splitmutex);

        // Doing the exciting bit now
        thisobj->splitBox(currentBoxCoords, particleData, targetData, IS_DYNAMIC_R);

        pthread_mutex_lock( &splitmutex); // prevent possible race condition to the while loop
    } // no more boxes on this level waiting to be split

    pthread_mutex_unlock( &splitmutex);
    return NULL;
}

//  ------------:
//  NAME        : .
//  PURPOSE     : Non-threaded standard box split algorithm.
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:
void BoxTree::normalBoxSplit(FmmTreeCoordList &splitBoxes, PointsContainer &particleData, PointsContainer &targetData, bool IS_DYNAMIC_R)
{
    for (FmmTreeCoordListIterator i = splitBoxes.begin(); i != splitBoxes.end(); i++)
    {
        splitBox(*i, particleData, targetData, IS_DYNAMIC_R);
    } // no more boxes on this level waiting to be split
}

//  ------------:
//  NAME        : .
//  PURPOSE     : Single box split algorithm.
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:
void BoxTree::splitBox(FmmTreeCoord boxCoords, PointsContainer &particleData, PointsContainer &targetData, bool IS_DYNAMIC_R)
{
    FmmBoxPointer currentBox = this->getBox(boxCoords);
    // If box has more vorts than 4*maxparticles then we will definately be splitting atleast one box...
    // but probably a waste of calculations as we can just do 3
    // extra times and will find same results later anyways

    // Particles that need to be sent to the children

    //	Now Cylce through all particles listed in this box, create children and update tree as neccesary and assign particles, keepng track of futuresplit list
    ComplexDouble particlePos, targetPos;
    DIRECTION_TYPE childDir;
    FmmTreeCoord childCoords;
    FmmBoxPointer childBox;
    size_t particleIndex, targetIndex;
    IntList * particles;
    IntList * targets;

    double newrsize, targetDist;

    double maxZdist[]={0,0,0,0};
   	particles = currentBox->getParticleList();

    ComplexDouble parentPos = currentBox->getPos();

    // Now cycle through particles, determining which box they fall under.
    for(IntListIterator i = particles->begin(); i != particles->end() ; i++)
    {
        particleIndex = *i;
        particlePos = particleData.getPos(particleIndex);

        // determine which child the particle would reside in
        childDir = calcChildDirection(parentPos, particlePos);

		if (!(FmmBoxClass::boxExists(currentBox->getChild(childDir))))
        {
            // this child didnt exist, so must create it before assigning particles to it
            createChild(currentBox, childDir);
        }

        childCoords = currentBox->getChild(childDir);   // resetting this variable to the valid particleData
        childBox = this->getBox(childCoords);
        // add particle to the child
        childBox->addParticle(particleIndex);
    } // allocated all particles now...technically they didnt need to be saved at any higher levels.

    // Now cycle through targets, determining which box they fall under.

    targets = currentBox->getTargetList();
    for(IntListIterator i = targets->begin(); i != targets->end() ; i++)
    {
        targetIndex = *i;
        targetPos = targetData.getPos(targetIndex);

        // determine which child the particle would reside in
        childDir = calcChildDirection(parentPos, targetPos);

		if (!(FmmBoxClass::boxExists(currentBox->getChild(childDir))))
        {
            // this child didnt exist, so must create it before assigning particles to it
            createChild(currentBox, childDir);
        }

        childCoords = currentBox->getChild(childDir);   // resetting this variable to the valid particleData
        childBox = this->getBox(childCoords);
        // add particle to the child
        childBox->addTarget(targetIndex);

        if (IS_DYNAMIC_R)
        {
            //determine if its further than the max distance!
            targetDist = abs(targetPos - childBox->getPos());
            if (targetDist > maxZdist[childDir])
            {
                maxZdist[childDir] = targetDist;
            }
        }
    } // allocated all particles now...technically they didnt need to be saved at any higher levels.
    //In an attempt to save memory this parents particles can now be removed!
    bool DELETE_AS_GO = true;
    if (DELETE_AS_GO)
    {
        //currentBox->setNumOfParticles(0);
        particles->clear();
        //delete particles;
        //currentBox->setNumOfTargets(0);
        targets->clear();
        //delete targets;
        //cout << "line 269 in treebuild changed";
    }
    // Fixing up the box rsizes according to the particle extremities
    if (IS_DYNAMIC_R)
    {
        for(size_t i = 0; i< MAX_CHILDREN; i++)
        {
            childCoords = currentBox->getChild(i);   // returns coords for the fmmtree
            if (FmmBoxClass::boxExists(childCoords))
            {
                // Doing this only on living children
                childBox = this->getBox(childCoords);
                // set rsize as max distance to particte
                newrsize=maxZdist[i];
                childBox->setRsize(newrsize);
            }
        }
    }
}


void * threadedPostSplit(void * inarguments)
{
    FmmTreeCoord currentBoxCoords;

    PostSplitArguments * arguments = static_cast<PostSplitArguments *>(inarguments);

	FmmTreeCoordList &splitBoxes = *arguments->splitBoxes;
	FmmTreeCoordList &futureSplitBoxes = *arguments->futureSplitBoxes;
   	size_t targetsPerBox = arguments->targetsPerBox;
   	size_t MIN_TREE_LEVEL = arguments->MIN_TREE_LEVEL;
	BoxTree * thisobj = arguments->thisobj;

    // prevent possible race condition to the while loop
    pthread_mutex_lock( &splitmutex);
    while (splitBoxIterator != splitBoxes.end())
    {
        // lets check us out a box in a thread friendly manner
        currentBoxCoords = *splitBoxIterator;
        splitBoxIterator++;
        pthread_mutex_unlock( &splitmutex);

        // Doing the exciting bit now
        thisobj->postSplit(currentBoxCoords, futureSplitBoxes, targetsPerBox, MIN_TREE_LEVEL);

        pthread_mutex_lock( &splitmutex); // prevent possible race condition to the while loop
    } // no more boxes on this level waiting to be split

    pthread_mutex_unlock( &splitmutex);
    return NULL;
}

//  ------------:
//  NAME        : .
//  PURPOSE     : Non-threaded standard box split algorithm.
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:
void BoxTree::normalPostSplit(FmmTreeCoordList &splitBoxes, FmmTreeCoordList &futureSplitBoxes, size_t targetsPerBox, size_t MIN_TREE_LEVEL)
{
    for (FmmTreeCoordListIterator i = splitBoxes.begin(); i != splitBoxes.end(); i++)
    {
        postSplit(*i, futureSplitBoxes, targetsPerBox, MIN_TREE_LEVEL);
    } // no more boxes on this level waiting to be split
}



//  ------------:
//  NAME        : .
//  PURPOSE     : Single post split algorithm.
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:
void BoxTree::postSplit(FmmTreeCoord boxCoords, FmmTreeCoordList &futureSplitBoxes, size_t targetsPerBox, size_t MIN_TREE_LEVEL)
{
    FmmBoxPointer currentBox = this->getBox(boxCoords);
    for(size_t i = 0; i< MAX_CHILDREN; i++)
    {
        FmmTreeCoord childCoords = currentBox->getChild(i);
        // returns coords for the fmmtree
        if (FmmBoxClass::boxExists(childCoords))
        {
            // Doing this only on living children
            FmmBoxPointer childBox = this->getBox(childCoords);
            setRelations(childBox, i, currentBox);

            // While we are at it, may aswell check if it needs to split!
            // Dont need to push_back particles if already at the maximimum level and the big loop isnt going to repeat!
            // Also if this box is less than the minimum tree level then all boxes will be
            bool treeIsTooShallow = (childBox->getLevel() < MIN_TREE_LEVEL);
            bool isTooManyTargets = (childBox->getNumTargets() > targetsPerBox);
            bool isTooManySources = (childBox->getNumParticles() > targetsPerBox);
            bool boxIsSplitting = (treeIsTooShallow || isTooManyTargets);

            if (!boxIsSplitting && isTooManySources)
            {
                // Putting in a little extra effort to see if we need to split due to neighbours
                size_t ineighb = 0;
                // Finding if neighbours should split
                FmmBoxPointer neighbourbox;
                FmmTreeCoord neighbourcoords;

                while (!boxIsSplitting && (ineighb < (MAX_NEIGHBOURS-1)))
                {
                    // Will stop searching once we either split or run out of neighbours
                    neighbourcoords = childBox->getNeighbour(ineighb);
                    if (FmmBoxClass::boxExists(neighbourcoords))
                    {
                        neighbourbox = this->getBox(neighbourcoords);
                        if (neighbourbox->getNumTargets() > targetsPerBox)
                        {
                            // Neighbour exists and has too many targets so it will be splitting
                            // This childbox has not enough targets to split, but has too many particles to allow it to become a phantom box.
                            // Dont think that i care if neighbour box splits due to its other neighbour interactions as it is outside this childbox's field of vision so it wont increase n^2 calcs at all.
                            boxIsSplitting = true;
                        }
                    }
                    ineighb++;
                }
            }


            if (boxIsSplitting)
            {
                //	Now if any child has too many vortexIndex add it to be split next time! In a threadsafe manner
                pthread_mutex_lock( &futureSplitMutex);
                futureSplitBoxes.push_back(childCoords);
                pthread_mutex_unlock( &futureSplitMutex);
            }
        }
    }
}




//  ------------:
//  NAME        : .
//  PURPOSE     : Sets up neighbours for all children.
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : Relies on the relations of this box to have been done at previous level
//  ------------:
void BoxTree::setRelations(FmmBoxPointer childBox, size_t childNum, FmmBoxPointer parentBox)
{
    FmmTreeCoord uncleCoords,neighbourCoord;
    FmmBoxPointer uncleBox;

    for(size_t ineighb = 0; ineighb < MAX_NEIGHBOURS; ineighb++)
    {
        uncleCoords = parentBox->getNeighbour( PARENT_TEMPLATE[childNum][ineighb] -1 ); // need c fix as template is defined in matlab indices
        if (FmmBoxClass::boxExists(uncleCoords))
        {
            uncleBox = this->getBox(uncleCoords);
            neighbourCoord = uncleBox->getChild( CHILD_TEMPLATE[childNum][ineighb] -1); // need c fix as template is defined in matlab indices

            childBox->setNeighbour(ineighb, neighbourCoord);
        }
    }
}

//  ------------:
//  NAME        : .
//  PURPOSE     : Uses number to add particles to the largest lev 0 box
//  IMPORTS     : .
//  EXPORTS     : .
//  PRE-CONDS   : .
//  POST-CONDS  : .
//  NOTES       : .
//  ------------:
FmmTreeCoord BoxTree::initializeFirstBox(PointsContainer &particleData, PointsContainer &targetData, const size_t inpValue)
{
    //	defing the size and location of the largest mesh.
    FmmTreeCoord god, boxCoords;
    double size;
    ComplexDouble centre;

    god.level = 0;
    god.index = 1;

    PointsContainer::getDomain(particleData, targetData, size, centre);

    // creating first box at level 0
    FmmBoxPointer firstbox = new FmmBoxClass(inpValue, god, centre, size);

    // adding it to the tree at level 0;
    this->initializeLevel(0,1);   //level 0 with 1 element
    boxCoords = addBox(0, firstbox);
    firstbox->setCoords(boxCoords);

    //need to seed itself in its neighbour relations
    firstbox->setNeighbour(SELF,boxCoords);

    // allocate all particles to it, as all are encapsulated at this level
    for(size_t i=0; i < particleData.getNumParticles(); i++)
    {
        firstbox->addParticle(i);
    }
    // allocate all targets to it, as all are encapsulated at this level
    for(size_t i=0; i < targetData.getNumParticles(); i++)
    {
        firstbox->addTarget(i);
    }

    return boxCoords;
}

void BoxTree::createChild(FmmBoxPointer parentBox, DIRECTION_TYPE childDir)
{
    // need to know new location of box centre first!
    ComplexDouble childPos = calcChildCentre(childDir, parentBox->getSize(), parentBox->getPos());

    FmmBoxPointer newbox = new FmmBoxClass(parentBox->getPvalue(), parentBox->getCoords(), childPos, ( parentBox->getSize() / 2 ));

    // adding it to the tree at level parent+1;
    pthread_mutex_lock( &addBoxMutex);
    FmmTreeCoord boxCoords = addBox(parentBox->getLevel()+1, newbox);
    pthread_mutex_unlock( &addBoxMutex);

    newbox->setCoords(boxCoords);

    // lastly must add child to the parents child field!
    parentBox->setChild(childDir, boxCoords);
}

DIRECTION_TYPE BoxTree::calcChildDirection(ComplexDouble parentPos, ComplexDouble particlePos)
{
    // //  simple function to determine the direction of a point relative to a
    // given centre
    DIRECTION_TYPE direction;

    if (real(particlePos) <= real(parentPos))
    {
        if (imag(particlePos) >= imag(parentPos))
        {
            direction = NW;
        }
        else
        {
            direction = SW;
        }
    }
    else
    {
        if (imag(particlePos) >= imag(parentPos))
        {
            direction = NE;
        }
        else
        {
            direction = SE;
        }
    }
    return direction;
}


ComplexDouble BoxTree::calcChildCentre(DIRECTION_TYPE boxDir, double parentSize, ComplexDouble parentCentre)
{
    double x,y,distToChild;

    distToChild = parentSize/4;

    if ( ( boxDir == NE ) || ( boxDir == NW ) )
    {
        y = imag(parentCentre) + distToChild;
    }
    else
    {
        y = imag(parentCentre) - distToChild;
    }

    if ( ( boxDir == NW ) || ( boxDir == SW ) )
    {
        x = real(parentCentre) - distToChild;
    }
    else
    {
        x = real(parentCentre) + distToChild;
    }

    return ComplexDouble(x,y);
}
