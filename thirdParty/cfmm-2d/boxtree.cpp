//==============|    
//  NAME        : fmmtreeclass.cpp
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 20.10.2008
//  DESCRIPTION : seperate implementaiton file for the FmmTree class, excluding the build function. 
//  NOTES       : Bug in maximum number of boxes limited to 256!
//  TODO        : 256 box max bug!! Heirachy functions for visualization!
//==============|
#include "boxtree.hpp"

#include <sstream>
#include <iostream>


BoxTree::BoxTree()
{
    fmmTree.clear();
    return;
}

BoxTree::~BoxTree()
{
    for(size_t i=0; i < fmmTree.size(); i++)
    {
        //cycle through each level, deleting the contents the levels point to.    
        for(size_t j=0; j < fmmTree[i]->size(); j++)
        {
            delete fmmTree[i]->operator[](j);
        }
        
        //then delete the level
        delete fmmTree[i];
    }
    return;
}








    
void BoxTree::initializeLevel(size_t level, size_t size)
{
    if (level < fmmTree.size())
    {
        cerr << "BoxTree: ERROR 1: oops tried to initialize the same level more than once" << endl;
    }      
    else
    {
        FmmBoxPointerArray* tmp = new FmmBoxPointerArray();
        tmp->reserve(size);
        fmmTree.push_back(tmp); 
    }
    return;
}

FmmTreeCoord BoxTree::addBox(size_t level, FmmBoxPointer box)
{
    FmmTreeCoord coords;
    size_t index;
    
    fmmTree[level]->push_back(box);
    index = fmmTree[level]->size();
    
    coords.level = level;
    coords.index = index;
    return coords;
}

size_t BoxTree::boxSizeToTreeLevel(double boxWidth, PointsContainer &particleData, PointsContainer &targetData )
{
    double domainsize;
    ComplexDouble centre;

    PointsContainer::getDomain(particleData, targetData, domainsize, centre);
    
    size_t treelevel = 0;
    if (domainsize > boxWidth)
    {
        treelevel = static_cast<size_t>(log(domainsize/boxWidth) / log(2));
    }
    return treelevel;    
}





const string BoxTree::toString(size_t lev, size_t index)
{
    
    ostringstream ss;

    FmmTreeCoord c;
    c.level = lev;
    c.index = index;
    FmmBoxPointer bx = getBox(c);
    
    ss << bx->toString() << bx->aktoString() << bx->bltoString()  << endl;
    return ss.str();

}
    
const string BoxTree::toString(size_t level)
{
    ostringstream ss;
    FmmBoxPointerArray* l = getLevel(level);
    

    ss << "Level " << level << " info----------" << endl;
    ss << "Address        = " << l << endl;
    ss << "numOfboxes     = " << l->size() << endl;
       
    if (l->size() <= 8)
    {
        ss << "LISTING BOXES: " << endl;
        for (size_t i =1; i <= l->size(); i++)
        {
             ss << toString(level,i) << endl;
        }
    }
    
	ss << endl;
	
    return ss.str();
}

const string BoxTree::toString()
{
    ostringstream ss;

    ss << "----------Tree Contains:----------" << endl;
    
    for (size_t i =0; i < fmmTree.size(); i++)
	{
	     ss << "^^^^^^^^^^^^^^^^^^" << endl << toString(i) << "vvvvvvvvvvvvvvvvvv" << endl;
    }
    return ss.str();
}

const string BoxTree::levelSizesToString()
{
    ostringstream ss;
    for (size_t i =0; i < fmmTree.size(); i++)
	{
        ss << " " << this->getNumBoxes(i) << ",";
        FmmTreeCoord master;
        master.level = i;
        master.index = 1;
        ss << endl <<  " LOWESTparts " << this->getBox(master)->getNumParticles();
        ss << " LOWESTparts_real " << this->getBox(master)->getParticleList()->size();
        ss << " LOWESTtargs " << this->getBox(master)->getNumTargets();
        ss << " LOWESTtargs_real " << this->getBox(master)->getTargetList()->size()<< endl;
    }
    return ss.str();
}


//~ 
//~ const string BoxTree::heirachyToString()
//~ {   
    //~ ostringstream ss;
//~ 
    //~ ss << "----------The FMM Family Tree:----------" << endl;
   //~ 
    //~ // reverse the tree starting from box (0,1)
    //~ 
    //~ FmmTreeCoord master;
    //~ master.level = 0;
    //~ master.index = 1;
    //~ 
    //~ FmmBoxPointer box = getBox(master);
    //~ ss << "starting at box coords " << box->getCoords() << endl;
    //~ 
    //~ ss << "(" << recursivePrint(box);       
    //~ 
    //~ return ss.str();
//~ }
//~ 
const string BoxTree::recursivePrint(FmmTreeCoord boxcoords)
{
    ostringstream ss; 
    //ostringstream str;
    //str << "str";
    
    //FmmTreeCoord coords = box->getCoords();
    FmmBoxPointer box = getBox(boxcoords);

    //str << boxcoords << " --> ";
    //~ str << box->getCoords() << ") --> (";
    
    //ss << ;
        
    //ss   
    if (!(box->isLowestBox()))
    {   
        FmmTreeCoord tmpboxcoord;
        //FmmBoxPointer newbox;
        for(size_t i=0; i < MAX_CHILDREN; i++)
        {
            tmpboxcoord = box->getChild(i);
            if (FmmBoxClass::boxExists(tmpboxcoord))
            {
                //newbox = getBox(tmpboxcoord);
                //ss << recursivePrint(newbox);
                ss << boxcoords << " --> " << recursivePrint(tmpboxcoord) << endl;
            }
        }
    }
    else
    {
        ss << endl;
    }
    //ss << str.str();

    return ss.str();    
}



ostream & operator<<(ostream &os, BoxTree &tree)
{
    os << tree.toString();
    return os;
}


//cvt pseudocode:
/*
 * 
 * create string
 * 
 * for level=0 to maxlevel
 * for box = 0 to maxbax
 *      // check for particles()
 *      if islowestbox
 *          get particles
 *          for i = 0 to numparticles
 *              checkparticlep[i] += 1
 *          
 *              if Particledistance_to_centre > boxSize 
 *                  ERROR particle in box not correct
 *          
 *       ifend
 *      // not lowest box so need to check children instead of particles
 *      for child = 1:4
 *          case 1,2,3,4) if childx,y in right quadrant and < parentsSize   
 *          else, children are out of whack!!
 *      end
 *endfor
 * endfor
 * 
 * for checkparticle = 0 to numparts
 * if cp[i] > 1 or == 0, ERROR missing box/counted twice 
 * 
 * 
 * 
 */

//~ const string BoxTree::confirmValidTree(PointsContainer &data)
//~ {
    //~ 
    //~ ostringstream ss;
//~ 
    //~ FmmBoxPointerArray* l;
    //~ FmmBoxPointer box, cbox;
    //~ IntList * plist;
    //~ double maxrdist, relx, rely;
    //~ ComplexDouble boxPos, rel, particlePos, vectToBox, cpos;
    //~ DIRECTION_TYPE dir;
    //~ FmmTreeCoord ne, sw, se, nw;
//~ 
    //~ size_t num=data.getNumParticles();
    //~ size_t checkparticles[num];
    //~ for(size_t i=0; i<num;i++)
    //~ {
       //~ checkparticles[i]=0;
    //~ }
//~ 
 //~ * for level=0 to maxlevel
    //~ for (size_t level =0; level < fmmTree.size(); level++)
	//~ {
	    //~ l= getLevel(level);
           //~ 
        //~ for (size_t i =0; i < l->size(); i++)
        //~ {
        //~ * for box = 0 to maxbax
              //~ box = l->at(i);
              //~ boxPos = box->getPos();
            //~ maxrdist = sqrt(2) * (box->getSize()/4);
  //~ *      // check for particles()
            //~ *      if islowestbox
 //~ *      
            //~ if (box->isLowestBox())
            //~ {
                //~ plist = box->getParticleList(); //   get particles
                //~ for (IntListIterator i =plist->begin(); i!= plist->end();i++)
                //~ {
                    //~ //for i = 0 to numparticles
                    //~ checkparticles[*i] += 1;    // count particle as being present
                    //~ //ss << "from lowest box" << box->getCoords() << ": incrementing Particle #" << *i << endl;
                                     //~ 
                    //~ particlePos = data.getPos(*i);
                    //~ vectToBox = particlePos - boxPos;
                    //~ // if Particledistance_to_centre > boxSize 
                    //~ if (abs(vectToBox) > box->getSize()) // !!!!!!!!!!!!! technically should be RSIZE!!!!!!
                    //~ {
                        //~ ss << "ERROR1 in box" << box->getCoords() << ": Particle #" << *i << " not enclosed within box!" << endl;
                    //~ }
                //~ }
            //~ }
            //~ else
            //~ {
                //~ // not lowest box so need to check children instead of particles
                //~ dir=NW;
                //~ nw = box->getChild(dir);    //for child = 1:4
                //~ if (FmmBoxClass::boxExists(nw))
                //~ {   
                    //~ cbox = getBox(nw);
                    //~ // is childbox in correct quadrant and within the parents size
                    //~ cpos = cbox->getPos();
                    //~ rel = cpos - boxPos;
                    //~ relx= real(rel);
                    //~ rely= imag(rel);
                    //~ 
                    //~ if ((relx > 0) || (rely <0))
                    //~ {
                         //~ ss << "ERROR2 in box" << box->getCoords() << " child(" << dir << ") is in wrong quadrant!" << endl;
                    //~ }
                    //~ 
                    //~ if (abs(rel) > maxrdist)
                    //~ {
                         //~ ss << "ERROR3 in box" << box->getCoords() << " child(" << dir << ") centre is too far away!" << endl;
                    //~ }
                //~ }
                //~ 
                //~ dir=SW;
                //~ sw = box->getChild(dir);    //for child = 1:4
                //~ if (FmmBoxClass::boxExists(sw))
                //~ {
                    //~ cbox = getBox(sw);
                    //~ // is childbox in correct quadrant and within the parents size
                    //~ cpos = cbox->getPos();
                    //~ rel = cpos - boxPos;
                    //~ relx= real(rel);
                    //~ rely= imag(rel);
                    //~ 
                    //~ if ((relx > 0) || (rely > 0))
                    //~ {
                        //~ ss << "ERROR2 in box" << box->getCoords() << " child(" << dir << ") is in wrong quadrant!" << endl;
                    //~ }
                    //~ 
                    //~ if (abs(rel) > maxrdist)
                    //~ {
                         //~ ss << "ERROR3 in box" << box->getCoords() << " child(" << dir << ") centre is too far away!" << endl;
                    //~ }
//~ 
                //~ }
                //~ 
                //~ dir=SE;
                //~ se = box->getChild(dir);    //for child = 1:4
                //~ if (FmmBoxClass::boxExists(se))
                //~ {
                    //~ cbox = getBox(se);
                    //~ // is childbox in correct quadrant and within the parents size
                    //~ cpos = cbox->getPos();
                    //~ rel = cpos - boxPos;
                    //~ relx= real(rel);
                    //~ rely= imag(rel);
                    //~ 
                    //~ if ((relx < 0) || (rely >0))
                    //~ {
                         //~ ss << "ERROR2 in box" << box->getCoords() << " child(" << dir << ") is in wrong quadrant!" << endl;
                    //~ }
                    //~ 
                    //~ if (abs(rel) > maxrdist)
                    //~ {
                         //~ ss << "ERROR3 in box" << box->getCoords() << " child(" << dir << ") centre is too far away!" << endl;
                    //~ }
                //~ }
                //~ 
                //~ dir=NE;
                //~ ne = box->getChild(dir);    //for child = 1:4
                //~ if (FmmBoxClass::boxExists(ne))
                //~ {
                    //~ cbox = getBox(ne);
                    //~ // is childbox in correct quadrant and within the parents size
                    //~ cpos = cbox->getPos();
                    //~ rel = cpos - boxPos;
                    //~ relx= real(rel);
                    //~ rely= imag(rel);
                    //~ 
                    //~ if ((relx < 0) || (rely <0))
                    //~ {
                         //~ ss << "ERROR2 in box" << box->getCoords() << " child(" << dir << ") is in wrong quadrant!" << endl;
                    //~ }
                    //~ 
                    //~ if (abs(rel) > maxrdist)
                    //~ {
                         //~ ss << "ERROR3 in box" << box->getCoords() << " child(" << dir << ") centre is too far away!" << endl;
                    //~ }
//~ 
                //~ }
            //~ }
        //~ }
    //~ }
//~ 
    //~ //lastly checking each particle was only counted once!
    //~ for (size_t i=0; i< data.getNumParticles(); i++)
    //~ {
        //~ if (checkparticles[i] != 1)
        //~ {
            //~ cerr << "ERROR4: Particle " << i << " count is equal to " << checkparticles[i] << endl;
        //~ }
    //~ }   
        //~ 
    //~ return ss.str();    
//~ }


//~ void interactivePrint()
//~ {
    //~ ostringstream ss;
    //~ //istringstream index;
    //~ 
    //~ size_t index, level;    
    //~ 
    //~ index=1;
    //~ cout << "Enter level, then index ( <= 0 to exit)" << endl;
    //~ cin >> level;
    //~ cin >> index;
    //~ 
    //~ while (index > 0) 
    //~ {
        //~ cout << "Enter level, then index ( <= 0 to exit)" << endl;
        //~ cin >> level;
        //~ cin >> index;  
        //~ coord.level=level;
        //~ coord.index=index;
        //~ 
        //~ parentsToString(coord)
        //~ 
                //~ 
    //~ }    
        //~ 
//~ }
