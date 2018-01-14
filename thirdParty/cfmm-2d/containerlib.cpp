//==============|    
//  NAME        : containerlib.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 07.10.2008
//  DESCRIPTION : Holds functions for data containers.
//  NOTES       : Relies on valid input data to ensure no seg faults and out of bound errors! 
//  TODO        : 
//==============|    
#include "containerlib.hpp"
//============================================================================
//=============== Constructors ===============================
//============================================================================
//~ PointsContainer::PointsContainer()
//~ {
//~ 
//~ }

PointsContainer::PointsContainer(double * inX, double *inY, double * inU, double *inV, size_t inNumParticles)
{
    if (inNumParticles > 0)
    {
        x = inX;
        y = inY;
        u = inU;
        v = inV;
        numParticles = inNumParticles;
        if ((x != NULL) && (y != NULL))
        {
            double minXa, maxXa, minYa, maxYa;
            calcExtremities(minXa, maxXa, minYa, maxYa, numParticles);
            setExtremities(minXa, maxXa, minYa, maxYa);
        }
    }
}

void PointsContainer::setPos(size_t partnum, ComplexDouble inPos)
{ 
    if (x == NULL || y == NULL)
    {
        cerr << "!!!x and y fields are null..." << endl;
    }
    else
    {
        x[partnum]= inPos.real();
        y[partnum]= inPos.imag();
    }
}

void PointsContainer::setVel(size_t partnum, ComplexDouble inVel)
{ 
    if (u == NULL || v == NULL)
    {
        cerr << "!!!u and v fields are null...meaning this isnt a targetclass!!" << endl;
    }
    else
    {
        u[partnum]= inVel.real();
        v[partnum]= inVel.imag();
    }
}


void PointsContainer::calcExtremities(double &minX, double &maxX, double &minY, double &maxY, size_t numParticles)
{
    // need to find max x, minx, maxy and minY       
    // easiest way is a forloop.
    double x, y;
    
    //need to seed it at a real point rather than 0,0! otherwise screws with boxcentre
    minX = getPos(0).real();
    maxX = minX;
    minY = getPos(0).imag();
    maxY = minY;    
    
    for (size_t i= 0; i < numParticles; i++)
    {   
        x = getPos(i).real();
        y = getPos(i).imag();
        
        if ( x < minX )
        {
            minX = x;
        }
        else if ( x > maxX )
        {
            maxX = x;
        }
        
        if ( y < minY )
        {
            minY = y;
        }
        else if ( y > maxY )
        {
            maxY = y;
        }
    } // found extremities of data
}

void PointsContainer::setExtremities(double minXa, double maxXa, double minYa, double maxYa)
{
    minX=minXa;
    maxX=maxXa;
    minY=minYa;
    maxY=maxYa;
}

void PointsContainer::getExtremities(double &outminX, double &outmaxX, double &outminY, double &outmaxY )
{
    outminX = minX;
    outmaxX = maxX;
    outminY = minY;
    outmaxY = maxY; 
}

double PointsContainer::getDomainSize()
{
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    double size;
    
    if (rangeX > rangeY)
    {
       size = rangeX;
    }
    else
    {
       size = rangeY;
    }
    return size;
}

ComplexDouble PointsContainer::getDomainCentre()
{
     return ComplexDouble( (maxX + minX)/2 , (maxY + minY)/2 );
}

void PointsContainer::getDomain(PointsContainer &particleData, PointsContainer &targetData, double &size, ComplexDouble &centre)
{
    double maxX, maxY,minX, minY, minXa, maxXa, minYa, maxYa ;
    
    particleData.getExtremities(minXa, maxXa, minYa, maxYa );
    targetData.getExtremities(minX, maxX, minY, maxY );
    
    // already defaulting to targets being the extremities
    if (minXa < minX)
    {
        minX = minXa;
    }    
    if (minYa < minY)
    {
        minY = minYa;
    }
    if (maxXa > maxX)
    {
        maxX = maxXa;
    }    
    if (maxYa > maxY)
    {
        maxY = maxYa;
    }
    
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    
    if (rangeX > rangeY)
    {
       size = rangeX;
    }
    else
    {
       size = rangeY;
    }
    
    centre = ComplexDouble( (maxX + minX)/2 , (maxY + minY)/2 );
    
}







//===========================================================================
//=============== ParticleContainer =================================        
//===========================================================================
//~ ParticleContainer::ParticleContainer()
//~ {
    //~ strength = NULL;
    //~ coreSize = NULL;
//~ }

ParticleContainer::ParticleContainer(double * inX, double *inY, double * inStrength, double * inCoreSizeSqrd, size_t inNumParticles)
:
PointsContainer(inX, inY, NULL, NULL, inNumParticles)
{
    // NOTE passing null for the u and v fields...as if we have a particle class we arent targetting it!
    if (inNumParticles > 0)
    {
        strength = inStrength;
        coreSizeSqrd = inCoreSizeSqrd;
    }
}

double ParticleContainer::getMaxSize()
{
    double maxsize,size; 
    
    //need to seed it at a real point rather than 0,0! otherwise its rooted
    maxsize = getCoreSizeSqrd(0);
    
    for (size_t i= 0; i < getNumParticles(); i++)
    {   
        size = getCoreSizeSqrd(i);
        if ( size > maxsize )
        {
            maxsize = size;
        }
    } // found extremities of data
    return maxsize;


}

//===========================================================================
//=============== ParticleMergeContainer =================================        
//===========================================================================
ParticleMergeContainer::ParticleMergeContainer(double * inX, double *inY, double * inStrength, double * inCoreSizeSqrd, size_t inNumParticles, double in_mincoreratio, double in_maxcoreratio, double in_radiusratio, double in_maxcoresizeratio)
:
    ParticleContainer::ParticleContainer(inX, inY, inStrength, inCoreSizeSqrd, inNumParticles)
{
    mincoreratio = in_mincoreratio;
    maxcoreratio = in_maxcoreratio;
    radiusratio = in_radiusratio;
    maxcoresizeratio = in_maxcoresizeratio;
    mergeCandidate = new bool[inNumParticles];

    for (size_t i = 0; i < inNumParticles; i++)
    {
        mergeCandidate[i] = true;
    }
}

ParticleMergeContainer::~ParticleMergeContainer()
{
    delete mergeCandidate;
}

void ParticleMergeContainer::getMergeParameters(double &in_mincoreratio, double &in_maxcoreratio, double &in_radiusratio, double &in_maxcoresizeratio)
{
    in_mincoreratio = mincoreratio;
    in_maxcoreratio = maxcoreratio;
    in_radiusratio = radiusratio;
    in_maxcoresizeratio = maxcoresizeratio;
}


//===========================================================================
//=============== SheetContainer =================================        
//===========================================================================
SheetContainer::SheetContainer(double * inXl, double *inYl, double * inStrengthl, double * inXr, double *inYr, double * inStrengthr, size_t numSheets)
:
ParticleContainer(NULL, NULL, NULL, NULL, numSheets)
{
    // NOTE passing null for all fields in the fathers container v fields...as if we have a particle class we arent targetting it!
    if (numSheets > 0)
    {
        leftSide = ParticleContainer(inXl, inYl, inStrengthl, NULL, numSheets);
        rightSide = ParticleContainer(inXr, inYr, inStrengthr, NULL, numSheets);        
    }
    double minXa, maxXa, minYa, maxYa;
    calcExtremities(minXa, maxXa, minYa, maxYa, numSheets);
    setExtremities(minXa, maxXa, minYa, maxYa);
}

double SheetContainer::getMaxLength()
{
    double maxsize,size; 
    
    //need to seed it at a real point rather than 0,0! otherwise its rooted
    maxsize = getLength(0);
    
    for (size_t i= 0; i < getNumParticles(); i++)
    {   
        size = getLength(i);
        if ( size > maxsize )
        {
            maxsize = size;
        }
    } // found extremities of data
    return maxsize;
}









//============================================================================
//=============== Other ===============================
//============================================================================

const string PointsContainer::toString(size_t partnum)
{
    ostringstream ss;
    ss << " pos=" << (getPos(partnum));
    if (u != NULL && v != NULL)
    {
        ss << " vel=" << (getVel(partnum)) ;
    }

    return ss.str();
}
const string PointsContainer::toString()
{
    string str;
    for(size_t i=0; i < getNumParticles() ; i++)
    {
        str += toString(i);
        str += "\n";
    }
    return str;
}

ostream & operator<<(ostream &os, PointsContainer &data)
{
    os << data.toString();
    return os;
}

const string ParticleContainer::toString(size_t partnum)
{
    ostringstream ss;
        ss << " str=" << getStrength(partnum);
        if (coreSizeSqrd != NULL)
        {
            ss << " size=" << (getCoreSizeSqrd(partnum));
        }
        ss << PointsContainer::toString(partnum);

    return ss.str();
}
const string ParticleContainer::toString()
{
    string str;
    for(size_t i=0; i < getNumParticles() ; i++)
    {
        str += toString(i);
        str += "\n";
    }
    return str;
}

ostream & operator<<(ostream &os, ParticleContainer &data)
{
    os << data.toString();
    return os;
}





const string SheetContainer::toString(size_t partnum)
{
    ostringstream ss;
    ss << " Length = " << getLength(partnum);
    ss << " |L|" << leftSide.toString(partnum);
    ss << " |R|" << rightSide.toString(partnum);

    return ss.str();
}
const string SheetContainer::toString()
{
    string str;
    for(size_t i=0; i < getNumParticles() ; i++)
    {
        str += toString(i);
        str += "\n";
    }
    return str;
}

ostream & operator<<(ostream &os, SheetContainer &data)
{
    os << data.toString();
    return os;
}
