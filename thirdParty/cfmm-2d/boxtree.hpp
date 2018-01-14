//==============|
//  NAME        : boxtree.hpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 20.10.2008
//  DESCRIPTION : Header for Fmmtree deals with a multilevelled structure of the and array of levels, where a level is an array of boxes.
//  NOTES       : uses STL
//  TODO        : Heirachy function for viewing
//==============|
#ifndef PARTICLECONTAINERCLASS_H
#define PARTICLECONTAINERCLASS_H
#include "containerlib.hpp"
#endif
#ifndef FMMBOXCLASS_H
#define FMMBOXCLASS_H
#include "fmmboxclass.hpp"
#endif

#ifndef _LIBHOARD
typedef vector<FmmBoxPointer> FmmBoxPointerArray;
#endif

#ifdef _LIBHOARD
typedef vector<FmmBoxPointer, std::allocator<FmmBoxPointer> > FmmBoxPointerArray;
#endif


//~ typedef vector<FmmBoxPointer, std::allocator<FmmBoxPointer>> FmmBoxPointerArray;
//    global    NW=1;NE=3;SE=4;SW=2;    N=5;E=6;S=7;W=8;SELF=9;
const size_t PARENT_TEMPLATE[4][9]= {
    {1,     8,     5,     9 ,    5 ,    9 ,    9 ,    8  ,   9},
    {8    , 2,     9,     7,     9,     9,     7,     8,     9},
    {5,     9,     3,     6,     5,     6,     9,     9,     9},
    {9,     7,     6,     4,     9,     6,     7,     9,     9}
                                    };
const size_t CHILD_TEMPLATE[4][9]=  {
    {4,     4,     4,     4,     2,     3,     2,     3,     1},
    {3,     3,     3,     3,     1,     4,     1,     4,     2},
    {2,     2,     2,     2,     4,     1,     4,     1,     3},
    {1,     1,     1,     1,     3,     2,     3,     2,     4}
                                    };

const size_t BOX_OFFSET = 1;

class BoxTree;

void * threadedBoxSplit(void *);
void * threadedPostSplit(void *);


typedef struct
{
	FmmTreeCoordList * splitBoxes;
	PointsContainer * particleData;
	PointsContainer * targetData;
	bool IS_DYNAMIC_R;
	BoxTree * thisobj;
} SplitBoxArguments;

typedef struct
{
	FmmTreeCoordList * splitBoxes;
	FmmTreeCoordList * futureSplitBoxes;
	size_t targetsPerBox;
	size_t MIN_TREE_LEVEL;
	BoxTree * thisobj;
} PostSplitArguments;






class BoxTree
{
public:
    BoxTree();
    ~BoxTree();

    FmmBoxPointerArray * getLevel(size_t level);
    FmmBoxPointer getBox(FmmTreeCoord coords);
    size_t getMaxLevel();
    size_t getNumBoxes(size_t lev);

    void initializeLevel(size_t level, size_t size);
	FmmTreeCoord addBox(size_t level, FmmBoxPointer box);
    void buildAdaptive(PointsContainer &particleData, PointsContainer &targetData, size_t inpValue, size_t targetsPerBox, size_t minTreeLevel, size_t maxTreeLevel,size_t maxThreads);

    static size_t boxSizeToTreeLevel(double boxWidth, PointsContainer &particleData, PointsContainer &targetData );
    const string toString(size_t lev, size_t index);
    const string toString(size_t level);
    const string toString();
    //~ const string heirachyToString();
    const string recursivePrint(FmmTreeCoord boxcoords);

    const string confirmValidTree(PointsContainer &data);
    const string levelSizesToString();

    friend void * threadedBoxSplit(void *);
    friend void * threadedPostSplit(void *);

private:
    vector<FmmBoxPointerArray*> fmmTree;

    //functions
    void normalBoxSplit(FmmTreeCoordList &splitBoxes, PointsContainer &particleData, PointsContainer &targetData, bool IS_DYNAMIC_R);
    void splitBox(FmmTreeCoord boxCoords, PointsContainer &particleData, PointsContainer &targetData, bool IS_DYNAMIC_R);

    void normalPostSplit(FmmTreeCoordList &splitBoxes, FmmTreeCoordList &futureSplitBoxes, size_t targetsPerBox, size_t MIN_TREE_LEVEL);
    void postSplit(FmmTreeCoord boxCoords, FmmTreeCoordList &futureSplitBoxes, size_t targetsPerBox, size_t MIN_TREE_LEVEL);

    FmmTreeCoord initializeFirstBox(PointsContainer &particleData, PointsContainer &targetData, const size_t inpValue);
    //~ void sendToChildren(PointsContainer &data, PointsContainer &targetData, FmmBoxPointer currentBox, bool IS_DYNAMIC_R);
    void createChild(FmmBoxPointer parentBox, DIRECTION_TYPE childDir);
    void setChildToSplit(FmmBoxPointer currentBox, FmmTreeCoordList &futureSplitBoxes);
    void calcRange(PointsContainer &data, double &size, ComplexDouble &centre);
    DIRECTION_TYPE calcChildDirection(ComplexDouble parentPos, ComplexDouble particlePos);
    ComplexDouble calcChildCentre(DIRECTION_TYPE boxDir, double parentSize, ComplexDouble parentCentre);
    void setRelations(FmmBoxPointer childBox, size_t childNum, FmmBoxPointer parentBox);
};

extern ostream & operator<<(ostream &os, BoxTree &tree);
//        void buildOldschool();



    //inline int getMaxLevel() { return maxTreeLevel;}		// Returns
    //inline int getMinLevel() { return minLevel;}			// what level is the highest that the upward pass went to.
    //inline double getLevelSize(level) {return (mesh0Size/(2^(level-1))); } 	// boxSize = mesh0Size/(2^(level-1));
    //inline double getLevelRsize(level) {return ((M_SQRT2/2.0) * getBoxsize(level)); } 	// boxRSize = (M_SQRT2/2.0) * mesh0Size/(2^(level-1));



//inlines
inline FmmBoxPointer BoxTree::getBox(FmmTreeCoord coords)
{
    // returns pointer to a sub-vector, then gets the box pointer from that vector object
    return fmmTree[coords.level]->operator[](coords.index-BOX_OFFSET);
}

inline FmmBoxPointerArray * BoxTree::getLevel(size_t level)
{
    return fmmTree[level];
}


inline size_t BoxTree::getMaxLevel()
{
    return (fmmTree.size()-1);
}

inline size_t BoxTree::getNumBoxes(size_t level)
{
    return fmmTree[level]->size();
}
