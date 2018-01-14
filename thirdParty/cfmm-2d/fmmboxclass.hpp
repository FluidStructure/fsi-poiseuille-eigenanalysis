//==============|    
//  NAME        : fmmboxclass.hpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 20.10.2008
//  DESCRIPTION : Headers to fmmboxclass.cpp, a class for... 
//  NOTES       : 
//  TODO        : add friend function for << operator on coords!
//==============|    

#ifndef LIST_H
#define LIST_H
#include <list>
#endif
#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif
#ifndef SSTREAM_H
#define SSTREAM_H
#include <sstream>
#endif
#ifndef IOSTREAM_H
#define IOSTREAM_H
#include <iostream>
#endif
#ifndef COMPLEX_H
#define COMPLEX_H
#include <complex>
#endif

using namespace std;
class FmmBoxClass;
struct FmmTreeCoord { size_t level; size_t index; };
typedef complex<double> ComplexDouble;
typedef FmmBoxClass * FmmBoxPointer;

#ifndef _LIBHOARD
typedef list<size_t> IntList;
typedef list<FmmTreeCoord> FmmTreeCoordList;
typedef vector<FmmTreeCoord> FmmTreeCoordArray;
typedef vector<size_t> IntArray;
typedef vector<ComplexDouble> ComplexDoubleArray;
#endif

#ifdef _LIBHOARD
typedef list<size_t, std::allocator<size_t> > IntList;
typedef list<FmmTreeCoord,std::allocator<FmmTreeCoord> > FmmTreeCoordList;
typedef vector<FmmTreeCoord,std::allocator<FmmTreeCoord> > FmmTreeCoordArray;
typedef vector<size_t,std::allocator<size_t> > IntArray;
typedef vector<ComplexDouble,std::allocator<ComplexDouble> > ComplexDoubleArray;
#endif

typedef IntList::iterator IntListIterator;
typedef FmmTreeCoordList::iterator FmmTreeCoordListIterator;

static const size_t MAX_CHILDREN=4;
typedef enum {NW=0, NE=2, SE=3, SW=1, N=4, E=5, S=6, W=7, SELF=8} DIRECTION_TYPE;
static const size_t MAX_NEIGHBOURS = 9;
static const FmmTreeCoord NULL_COORDS = { 0, 0 };

extern ostream &operator<<(ostream &os, FmmTreeCoord const &coord);

class FmmBoxClass 
{
// Public API
public:

    // Con/Destructors
    FmmBoxClass(size_t inpValue, FmmTreeCoord inparent, ComplexDouble pos, double insize);
    ~FmmBoxClass();

    // Accessors
    // multipole stuff
    ComplexDouble getAkCoeff(size_t pVal); 
    ComplexDouble getBlCoeff(size_t pVal);
    
    // box
    size_t getNumParticles();
    size_t getNumTargets();
    ComplexDouble getPos();
    size_t getLevel();
    double getSize();
    size_t getPvalue();
    FmmTreeCoord getCoords();
    FmmTreeCoord getParent();
    FmmTreeCoord getChild(size_t dir);
    FmmTreeCoord getNeighbour(size_t dir);
    size_t getNumPhantomBoxes();
    IntList * getParticleList();
    IntList * getTargetList();
    FmmTreeCoordList * getPhantomBoxList();
    double getRsize();
 
    // Mutators    
    void setCoords(FmmTreeCoord coords);
    void setRsize(double rSize);

    void setAkCoeff(size_t pVal, ComplexDouble inAkCoeff); // need to add check if not null then theres memory leaking!
    void setBlCoeff(size_t pVal, ComplexDouble inBlCoeff); // need to add check if not null then theres memory leaking!
    void addPhantomBox(FmmTreeCoord box);
    void addParticle(size_t particleInd);
    void addTarget(size_t particleInd);
    
    void setChild(size_t childnum, FmmTreeCoord coords);
    void setPhantomBoxlist(FmmTreeCoordList * phantomBoxList);
    void setNeighbour(size_t neighbnum, FmmTreeCoord neighbourCoord);

    // Methods
    static bool boxExists(FmmTreeCoord box);
    bool particleExists(size_t particle);
    bool isWellSeperated(FmmBoxPointer otherBox);
    bool isLowestBox();
    bool holdsTargets(); // This is always true...never split a box unless there is targets in it??? no, not true, can split if its less than split mark. always true
    bool holdsParticles();
           
    const string toString();
    const string aktoString();
    const string bltoString();


      
// Private fields/functions
private:
	double size, rSize;
	FmmTreeCoord coords;
	ComplexDouble centerPos;
	FmmTreeCoord parent;
    FmmTreeCoord * childList;
    FmmTreeCoord neighbourList[MAX_NEIGHBOURS];
    
    ComplexDouble * aKCoeff;
    ComplexDouble * bLCoeff; 
    IntList particleList;
    IntList targetList;
    FmmTreeCoordList phantomBoxList;
    size_t pvalue, numOfParticles, numOfTargets;

};



//---inline functions 
//============================================================================
//=============== Accessors ===============================
//============================================================================
    
inline ComplexDouble FmmBoxClass::getPos() 
{
    return centerPos;
}

inline FmmTreeCoord FmmBoxClass::getCoords()
{
    return coords;
}

inline size_t FmmBoxClass::getPvalue()
{
    return (pvalue);
}

inline double FmmBoxClass::getSize()
{
    return size;
}

inline ComplexDouble FmmBoxClass::getAkCoeff(size_t pval)
{
    return aKCoeff[pval]; 
}

inline ComplexDouble FmmBoxClass::getBlCoeff(size_t pval)
{
    return bLCoeff[pval]; 
}

inline double FmmBoxClass::getRsize()
{
    return rSize;
}

inline size_t FmmBoxClass::getNumPhantomBoxes()
{
    return phantomBoxList.size();
}

inline FmmTreeCoord FmmBoxClass::getNeighbour(size_t neighbnum)
{
    return neighbourList[neighbnum];
}
    
inline FmmTreeCoord FmmBoxClass::getParent()
{
    return parent;
}






inline IntList * FmmBoxClass::getParticleList()
{
    return &particleList;
}

inline IntList * FmmBoxClass::getTargetList()
{
    return &targetList;
}

inline size_t FmmBoxClass::getNumParticles()
{
    return this->numOfParticles;
}

inline size_t FmmBoxClass::getNumTargets()
{
    return this->numOfTargets;
}

inline size_t FmmBoxClass::getLevel()
{
    return coords.level;
}

inline FmmTreeCoordList * FmmBoxClass::getPhantomBoxList()
{
    return &phantomBoxList;
}

//============================================================================
//=============== Mutators ===============================
//============================================================================

inline void FmmBoxClass::setCoords(FmmTreeCoord incoords)
{
    coords=incoords;
}

inline void FmmBoxClass::setRsize(double inrSize)
{
    rSize=inrSize;
}

    

inline void FmmBoxClass::setAkCoeff(size_t pVal, ComplexDouble inAkCoeff)
{
    aKCoeff[pVal] = inAkCoeff;
}

inline void FmmBoxClass::setBlCoeff(size_t pVal, ComplexDouble inBlCoeff) 
{
    bLCoeff[pVal] = inBlCoeff;
}

inline void FmmBoxClass::addPhantomBox(FmmTreeCoord box)
{
    phantomBoxList.push_back(box);
}

inline void FmmBoxClass::addParticle(size_t particleInd)
{   
    particleList.push_back(particleInd);
    numOfParticles++;
}

inline void FmmBoxClass::addTarget(size_t particleInd)
{   
    targetList.push_back(particleInd);
    numOfTargets++;
}

inline void FmmBoxClass::setNeighbour(size_t neighbnum, FmmTreeCoord neighb)
{
    neighbourList[neighbnum] = neighb;
}

//============================================================================
//=============== Other ===============================
//============================================================================    

inline bool FmmBoxClass::isLowestBox()
{
    return (childList == NULL);
}
inline bool FmmBoxClass::holdsTargets()
{
    return (getNumTargets() > 0);
}
inline bool FmmBoxClass::holdsParticles()
{
    return (getNumParticles() > 0);
}

inline bool FmmBoxClass::boxExists(FmmTreeCoord coords)
{
    return (coords.index >0);// && coords.level >= 0 );
}
