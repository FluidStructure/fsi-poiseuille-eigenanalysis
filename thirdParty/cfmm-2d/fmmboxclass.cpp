//==============|    
//  NAME        : fmmboxclass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 23.10.2008
//  DESCRIPTION : Class to handle a single FMM box.
//==============|    
#include "fmmboxclass.hpp"

//============================================================================
//=============== Constructors ===============================
//============================================================================
FmmBoxClass::FmmBoxClass(size_t inpValue, FmmTreeCoord inparent, ComplexDouble pos, double insize) 
{
    size = insize;
    rSize= sqrt(2.0) * size / 2.0;
    centerPos = pos;
    pvalue=inpValue;
    
    if (inpValue > 0)
    {
        aKCoeff = new ComplexDouble[inpValue+1];
        bLCoeff = new ComplexDouble[inpValue+1];
        ComplexDouble zero = ComplexDouble(0);
        for (size_t i=0; i<=inpValue;i++)
        {
            aKCoeff[i]=zero;
            bLCoeff[i]=zero;
        }
    }

    for (size_t i=0; i < MAX_NEIGHBOURS;i++)
    {
        neighbourList[i]=NULL_COORDS;
    }
	childList=NULL;
      
    particleList = IntList();
    targetList = IntList();
    numOfParticles = 0;
    numOfTargets = 0;
    phantomBoxList = FmmTreeCoordList();
    parent=inparent;
}

FmmBoxClass::~FmmBoxClass()
{
    if (pvalue>0)
    {
    delete [] aKCoeff;
    delete [] bLCoeff;
    }
    if (!(isLowestBox()))
    {
    delete [] childList;
    }
}

FmmTreeCoord FmmBoxClass::getChild(size_t dir)
{
    FmmTreeCoord coords;
    if (childList==NULL)
    {
        coords=NULL_COORDS;
    }
    else
    {
        coords = childList[dir]; // dir should be automatically converted to its int value
    }
    return coords;
}



void FmmBoxClass::setChild(size_t child, FmmTreeCoord coords) 
{
    if (childList==NULL)
    {
        childList=new FmmTreeCoord[MAX_CHILDREN];
        for (size_t i=0; i<MAX_CHILDREN;i++)
        {
            childList[i]=NULL_COORDS;
        }
    }
    childList[child] = coords;
}


void FmmBoxClass::setPhantomBoxlist(FmmTreeCoordList * inphantomBoxList)
{
    if (phantomBoxList.size() > 0)
    {
        cerr << "DELETING A PHANTOM BOXLIST WHICH != 0 size" << endl;
    }            
    phantomBoxList = FmmTreeCoordList(*inphantomBoxList);
}

bool FmmBoxClass::isWellSeperated(FmmBoxPointer boxtwo)
{
    bool wellsep;
    FmmTreeCoord compare,neighb;
    
    wellsep = true;
    compare = boxtwo->getCoords();
    for (size_t n=0;n < MAX_NEIGHBOURS; n++)
    {
        neighb = getNeighbour(n);
        if (neighb.index == compare.index)
        {
            wellsep = false;
        }
    }
    
    //~ (abs(getPos() - boxtwo->getPos()) > 3.0 * getSize());
    return wellsep;
}

// IO funs use for debugging
const string FmmBoxClass::aktoString()
{
    ostringstream ss;
    ss << "ak array equal to --------: " << endl;    
    for (size_t i=0; i < (getPvalue() + 1); i++)
	{   
        ss << i << " : " << getAkCoeff(i) << endl;
    }
    return ss.str();
}

const string FmmBoxClass::bltoString()
{
    ostringstream ss;
    ss << "bl array equal to --------: " << endl;    
    for (size_t i=0; i < getPvalue() + 1; i++)
	{
    	ss << i << " : " << getBlCoeff(i) << endl;
    }
    return ss.str();
}

const string FmmBoxClass::toString()
{   
    ostringstream ss;

    ss << "Displaying info for box stored at " << this << "----------" << endl;
    ss << "Tree Coords    = " << getCoords() << endl;
    ss << "BoxSize        = " << getSize() << endl;
    ss << "centerPos      = " << getPos() << endl;
    ss << "parent address = " << getParent() << endl;
		
	ss << "Child addresses= ";

	if (isLowestBox())
	{
	    ss << "None";
    }
    else
    {
        for (size_t i = 0; i <MAX_CHILDREN; i++)
        {
                ss << getChild(i) << " , ";
        }
    }
	ss << endl;
	
	ss << "Neigb addresses= ";
	for (size_t i = 0; i <MAX_NEIGHBOURS; i++)
	{
	    ss << getNeighbour(i) << " , " ;
    }
	ss << endl;
	
	
	IntList * particleList = getParticleList();
	ss << "Particle list  = ";
	if (particleList->size() <=10000 )
    {    
        for(IntListIterator i= particleList->begin(); i != particleList->end(); i++)
        {
            ss << *i << " , ";
        }    
        ss << endl;
    }
    else
    {
        ss << "TOO MANY" << endl;
    }
    
    
    ss << "target list  = ";
	IntList * targetList = getTargetList();
	if (targetList->size() <=10 )
    {
        for(IntListIterator i= targetList->begin(); i != targetList->end(); i++)
        {
            ss << *i << " , ";
        }    
        ss << endl;
    }
        else
    {
        ss << "TOO MANY" << endl;
    }
	
	
    FmmTreeCoordList * phantomBoxList = getPhantomBoxList();
    ss << "Phantom addresses= ";
    for(FmmTreeCoordListIterator i= phantomBoxList->begin(); i != phantomBoxList->end(); i++)
	{
	    ss << *i << " , ";
    }     
    ss << endl;
    return ss.str();
}





ostream & operator<<(ostream &os,const FmmTreeCoord &coord)
{
    os << "(" << coord.level << "," << coord.index << ")";
    return os;
}
