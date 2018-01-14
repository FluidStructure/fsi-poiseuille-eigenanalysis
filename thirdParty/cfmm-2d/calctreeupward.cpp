//==============|    
//  NAME        : calcupwardpass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 07.10.2008
//  DESCRIPTION : Upward pass of a fast multipole algorithm adapted from (Greengard, 1987) fast alg for coulombic particle systems
//  NOTES       : Is not a standalone file, requires fastmultipolemethod function to call it into action
//  TODO        : 
//==============|    

#include "treesolverlib.hpp"

pthread_mutex_t upwardLevelMutex = PTHREAD_MUTEX_INITIALIZER;
size_t upwardBoxCounter = 0; //used for the threads checking out their own box...must be reset every level 

//  ------------:
//  NAME:         calcAkCoeff
//  PURVelE     : Calculates the A_k values from a single point to its containing box centre.
//  IMPORTS     : fmmTree object
//  PRE-CONDS   : Requires a valid fmmTree to work on.
//  POST-CONDS  : outAkCoeff will be set for all boxes in the fmmtree
//  NOTES       : .
//  ------------:
void FastMultipoleMethod::calcTreeUpward(BoxTree &fmmTree, ParticleContainer &data, vector<IntArray*> &binarray, size_t maxThreads)
{    
    size_t THREADS=maxThreads;
    size_t MIN_TREE_LEVEL = 2; //fmm shouldnt go any higher than this!
    // Cycle through each level of the quad tree from bottom up and calculate the actual and translated multipole expansions due to the contained particles. Only going up to level 2, as can't use any expansions on level 0 or 1 as no boxes will be "well seperated".
    for (size_t ilevel = fmmTree.getMaxLevel(); ilevel >= MIN_TREE_LEVEL; ilevel--)
    {
        // Start at level 1 because multipole method doesnt work at level 1 (points no longer are well seperated)
        // otherwise the children wont have a bl set    
        upwardBoxCounter = 0;     // must be reset every level

        if (THREADS > 1)
	    {
            pthread_t threadID[THREADS];
            // sorting out input arguments
            UpwardArguments threadArgs; 
            
            threadArgs.level = ilevel;
            threadArgs.fmmTree = &fmmTree;
            threadArgs.particleData = &data;
            threadArgs.binarray = &binarray;
            threadArgs.thisobj = this;

            for(size_t i=0; i < THREADS; i++)
            {                
                pthread_create( &threadID[i], NULL, threadedLevelUpward, (void *) &threadArgs );
            }
            //make sure threads have joined before moving on
            for(size_t i=0; i < THREADS; i++)
            {
                pthread_join( threadID[i], NULL); 
            }
        }
        else
        {
            normalLevelUpward(ilevel,fmmTree,data,binarray);
        }        
        
    } // Upward pass completed. There wont be any unevaluated Ak coefficients
} 

void * threadedLevelUpward(void * inarguments) 
{
    FmmTreeCoord boxcoords;
    UpwardArguments * arguments = static_cast<UpwardArguments *>(inarguments);
  
    size_t ilevel = arguments->level;   
    ParticleContainer &particleData = *arguments->particleData;
	BoxTree & fmmTree= *arguments->fmmTree;
	vector<IntArray*> binarray = *arguments->binarray;
	FastMultipoleMethod * thisobj = arguments->thisobj;
		
    boxcoords.level = ilevel;

	size_t totalBoxes = fmmTree.getNumBoxes(ilevel);
    pthread_mutex_lock( &upwardLevelMutex);
    while (upwardBoxCounter < totalBoxes)
    {   
        // cycle through all boxes on this level
        // checkout this box so no other threads try and do the same one!
        upwardBoxCounter++;
        boxcoords.index = upwardBoxCounter;
        pthread_mutex_unlock( &upwardLevelMutex);
        //other threads are now able to continue the loop
        
        // Doing the exciting bit now
        thisobj->levelUpward(boxcoords, fmmTree, particleData, binarray);
        
        pthread_mutex_lock( &upwardLevelMutex);
    } // No more boxes on this level
    pthread_mutex_unlock( &upwardLevelMutex);

    return NULL;
}  
 
void FastMultipoleMethod::normalLevelUpward( size_t ilevel, BoxTree &fmmTree, ParticleContainer &particleData,vector<IntArray*> &binarray)
{
    FmmTreeCoord boxcoords;

    boxcoords.level=ilevel;
    // Cycling through every box on this level
    for (size_t ibox = 1; ibox <= fmmTree.getNumBoxes(ilevel); ibox++)
    {
        boxcoords.index=ibox;
        levelUpward( boxcoords, fmmTree, particleData, binarray);
     } // Level finished, multipole calculated for all boxes on this level
}
 
void FastMultipoleMethod::levelUpward( FmmTreeCoord boxcoords, BoxTree &fmmTree, ParticleContainer &particleData, vector<IntArray*> &binarray)
{
    FmmBoxPointer box = fmmTree.getBox(boxcoords);
    if (box->holdsParticles())
    {
        // Only do upward motion for boxes with particles in it                   
        if (box->isLowestBox())
        {
            // There is no children to translate a multipole from so it must then be a particle holding box which we will do an actual multipole expansion about its centre to get A_k. 
            calcMultipoleExpansion(box, particleData);                    // The A_k multipole expansions truncated to pvalue has been completed at this lowest box.    Will be performing only local translations for the remaining traversing of this branch.                  
        } // Box finished, A_k now valid.
        else
        {
            // This box has no actual particles, but its descendants do. All childrens A_k expansions were done at the previous level so only thing to do is translate the expansion from all 4 children to parent. 
            translateMultipoleExpansions(box, fmmTree, binarray);
        }
    } // Box finished, contains a multipole expansion.
}

 
 
  
//============================================================================
//=============== FOR MULTIPOLE CALCULATION ==================================        
//============================================================================

//  ------------:
//  NAME        : calcMultipoleExpansion
//  PURPOSE     : Calculates the multipole expansion of a whole box for all the
//      vortices contained in it.
//  IMPORTS     : box pointer
//                the particleContainer
//              pValue
//  PRE-CONDS   : box is a lowest level box needing a multipole expansion
//  POST-CONDS  : boxes ak field is set and valid.
//  ------------:                    
void FastMultipoleMethod::calcMultipoleExpansion(FmmBoxPointer box, ParticleContainer &data)
{
    ComplexDouble ak,zRel,boxPos,particlePos;
    double particleStr;
    size_t maxPval;
    IntList * particleList;
    
    maxPval = box->getPvalue();
    boxPos = box->getPos();
    // cycle through all particles in this box:
    particleList = box->getParticleList();
    
    for (IntListIterator i=particleList->begin(); i != particleList->end(); i++)
    {    
        // For each particle in box
        particleStr = data.getStrength(*i);
        particlePos = data.getPos(*i);
        zRel = particlePos - boxPos;
        
        ak = particleStr + box->getAkCoeff(0);
        box->setAkCoeff(0 , ak); //  A_0 is the particle strength
        // for the rest of the A_ values
        for(size_t k=1; k <= maxPval; k++)
        {    
            // from lemma 2.? in Greengard
            ak = -particleStr * pow(zRel,k) / (double)k;
            
            ak += box->getAkCoeff(k);
            box->setAkCoeff(k , ak);
        } // the aKcoeff of the multipole expansion calculated for this particle
    } // Expansion of all particles completed, hence Akcoeff now set for this box.                                          
    return;
}

//============================================================================
//=============== FOR EXISTING MULTIPOLE TRANSLATIONS ======================== 
//============================================================================

//  ------------:
//  NAME        : translateMultipoleExpansion
//  PURPOSE     : Translates the multipole at a single childs center onto its parents center.
//  IMPORTS     : fmmTree (needed for out of box/level access)
//  PRE-CONDS   : bring in a box whose children have valid multipole AK coeffs
//  POST-CONDS  : .
//  NOTES       : .
//  ------------:
void FastMultipoleMethod::translateMultipoleExpansions(FmmBoxPointer box, BoxTree &fmmTree, vector<IntArray*> &binarray)
{
    ComplexDouble ak,boxPos,childPos,childQ,childAk,zRel;
    FmmTreeCoord childCoord;    
    FmmBoxPointer childBox;
    size_t pval;
    
    pval = box->getPvalue();
    boxPos = box->getPos();
    
    for(size_t i = 0; i < MAX_CHILDREN; i++)
    {
        // Cycling through all children
        childCoord = box->getChild(i);                                
        if (FmmBoxClass::boxExists(childCoord))
        {    
            // This child actually had particles in it an hence got created/exists.
            childBox = fmmTree.getBox(childCoord);  // Real memory address of child box
            if (childBox->holdsParticles())
            {
                // Only do upward motion for children with particles in it    
                childPos = childBox->getPos();
                childQ = childBox->getAkCoeff(0);
                zRel = childPos - boxPos;
                
                ak = box->getAkCoeff(0) + childQ;                                
                box->setAkCoeff(0, ak);                                

                for(size_t l=1; l <= pval; l++)
                {
                    ak = -childQ * pow(zRel,l) / (double)l;
                    for(size_t k=1; k <= l; k++)
                    {    
                        //    from lemma 2.3 in Greengard
                        childAk = childBox->getAkCoeff(k);
                        ak += childAk * pow(zRel,(l-k)) * bincoeff(l-1,k-1, binarray);  
                    } // the aKcoeff of the multipole expansion calculated for this particle
                    
                    ak += box->getAkCoeff(l);       // basically incrementing ak
                    box->setAkCoeff(l, ak);
                }   
            }
        }
    }
    return;
}
