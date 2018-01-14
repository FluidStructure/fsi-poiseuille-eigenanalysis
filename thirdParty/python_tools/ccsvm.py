#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Contains ccsvm functions.

AUTHOR - Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
MODIFIED - 2009.12
"""
import pdb
import math
import random
from velocitylib import Elemat
from velocitymodule import lamb_vortex_merge as lamb_vortex_merge_cfmm
def lamb_vortex_split(vormat, split_inds, num_children, split_alpha, split_scheme='huang4', position_scheme='random', blsim_dict=None):
    """For CCSVM splitting.
    
    USAGE         vormat = matrix containing discrete vortices in form [vortX, vortY, vortStr, vortCoreSizeSqrd]
                  maxcoresqrdsize = maximum sigma^2 value before splitting
                  splitalpha = the accuracy of the splitting (number between 1 and 0) closer to 1 the more accurate the spliting algorithm. However new children cores will almost be as big as parents
                  num_children = number of children to spawn around the parent
                  splitScheme = the type of splitting method to use.
                                1 : scheme IV in Huang, 2005
    -----------
    
    """

    if blsim_dict is not None:
        from boundarytheorylib import calc_lam_bound_lay_vel
        near_wall_cuttoff = 0.0 # not really true...but too hard with a moving wall! not really, just pass it in.!
        #bl_sim=False, profile_core_size_buffer=None, 
        bndary_height=blsim_dict['bndary_height']
        orig_bndary_height = bndary_height
        find_height_dict = blsim_dict['find_height_dict']
        find_height_fun = find_height_dict['fun']
        wallvortoffset = find_height_dict['offset']
        wallvormat = find_height_dict['wallvormat']
        find_optimum_coresize = blsim_dict['find_optimum_coresize']
        error_ratio = blsim_dict['error_ratio']
        core_edge_detect_ratio = blsim_dict['core_edge_detect_ratio']
    #None, core_edge_detect_ratio=None, error_ratio=None, prev_vormat_u=None, prev_vormat_v=None, b 
            # How small discretisation of domain is for stepping through, dont need to much precision as its a semi-fudge anyways.
    
    # Find vortices that have grown too large.
    for i in split_inds:
        split_x = vormat['x'][i]
        split_y = vormat['y'][i]
        split_s = vormat['s'][i]
        split_c = vormat['c'][i]

        if (split_scheme == 'huang4'):
            # represents scheme IV in Huang, 2005
            # never need to delete vortices...they can merge if need be      
            # Determining locations of children
            childradius = math.sqrt((split_c ** 2) * 2.0 * (1.0 - split_alpha ** 2))

            # Determining strength
            childstrength = split_s / (2.0 * num_children)

            # Determining coresize
            childcoresize = split_alpha * split_c
            if blsim_dict is not None: 
                # first test that the first guess wasnt adequite...
                # worst case scenario is a core directly down...maybe can sense this better judging 
                # this is a little bit of an overkill, as it could be at best 45 degrees off this
                real_center = split_y - childradius
                wall_height = find_height_fun(wallvormat, split_x, wallvortoffset, relative_to_sublayer=True)
                if real_center < wall_height:
                    childcoresize = 0.0
                else:
                    childcoresize = find_optimum_coresize(split_x, real_center, error_ratio, core_edge_detect_ratio, bndary_height, wallvormat, wallvortoffset, childcoresize)

            # Fixing parents: Location stays same, strength is halfed, coresize
            # is set to same as its children
            vormat['s'][i] = split_s / 2.0
            vormat['c'][i] = childcoresize

            # Initialise the new children
            # Get the angle between each child
            childseperation = 2.0 * math.pi / num_children
            # Get a random starting angle
            if position_scheme == 'random':
                childangle = random.random() * 2.0 * math.pi
            else:
                raise ValueError('Requested split position scheme' + str(position_scheme) + 'not recognised...!')
                
            # Inserting children into vormat
            # Add the new vortex to the end of vormat...is safe ONLY because nothing is being deleted in front of it
            for i in xrange(num_children):
                # Get the angle to put the new vortices in
                ang = childangle + childseperation * i
                # Get the x and y pos using the "sin" and "cos" of the angle
                thischildx = split_x + math.cos(ang) * childradius
                thischildy = split_y + math.sin(ang) * childradius
                # Add the new vortex to the end of vormat...is safe ONLY because nothing is being deleted in front of it
                vormat.add_particle(thischildx, thischildy, childstrength, childcoresize)
        
        elif (split_scheme == 'rossi'):
            raise Exception("rossi split not implemented yet")
                        
            # Determine coresize first
            # Determining locations of children
            childcoresize = split_alpha * split_c
            childradius = childcoresize * math.sqrt(1.0 - split_alpha ** 2)

            # Determining strength
            childstrength = split_s / num_children

            # Fixing parents: Location stays same, strength is halfed, coresize
            # is set to same as its children
            vormat['s'][i] = split_s / 2.0
            vormat['c'][i] = childcoresize

            # Initialise the new children
            # Get the angle between each child
            childseperation = 2.0 * math.pi / num_children
            # Get a random starting angle
            if position_scheme == 'random':
                childangle = random.random() * 2.0 * math.pi
            else:
                raise ValueError('Requested split position scheme' + str(position_scheme) + 'not recognised...!')
                
            # Inserting children into vormat
            # Add the new vortex to the end of vormat...is safe ONLY because nothing is being deleted in front of it
            for i in xrange(num_children):
                # Get the angle to put the new vortices in
                ang = childangle + childseperation * i
                # Get the x and y pos using the "sin" and "cos" of the angle
                thischildx = split_x + math.cos(ang) * childradius
                thischildy = split_y + math.sin(ang) * childradius
                # Add the new vortex to the end of vormat...is safe ONLY because nothing is being deleted in front of it
                vormat.add_particle(thischildx, thischildy, childstrength, childcoresize)
              
        else:
            raise ValueError('Requested splitting scheme \"' + str(split_scheme) + '\" not recognised...!')


def lamb_vortex_merge(vormat, mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio, fmm=True, threads=4, prev_dict=None):
    """For CCSVM merging in C++f
    
    Suggested parameters 0.9, 1.1, 0.5, 2.0
    
    """
    # doesnt matter, algorithm not multilevelled...
    fmm_targets_per_box = 15
    
    # FMM likes lists only! 
    vort_x = vormat['x']
    vort_y = vormat['y']
    vort_str = vormat['s']
    vort_core_size = vormat['c']

    if not hasattr(vort_x, '__len__'):
        vort_x = [vort_x]
        vort_y = [vort_y]
        vort_str = [vort_str]
        vort_core_size = [vort_core_size]
   
    if prev_dict is None:
        prev_u = [0.0] * len(vormat['x'])
        prev_v = [0.0] * len(vormat['x'])
    else: 
        prev_u = prev_dict['u']
        prev_v = prev_dict['v']

    out_x = list()
    out_y = list()
    out_s = list()
    out_c = list()
    out_pu = list()
    out_pv = list()


    lamb_vortex_merge_cfmm(mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio, fmm, fmm_targets_per_box, threads, vort_x, vort_y, vort_str, vort_core_size, prev_u, prev_v, out_x, out_y, out_s, out_c, out_pu, out_pv)

    # now pass the new vormat data back
    vormat['x'] = out_x
    vormat['y'] = out_y
    vormat['s'] = out_s
    vormat['c'] = out_c
    
    if prev_dict is not None:
        # not really robust but it will have to do
        # replace 0.0 with nones
        for i in xrange(len(out_pu)):
            if out_pu[i] == 0.0 and out_pv[i] == 0.0:
                out_pu[i] = None
                out_pv[i] = None
        
        prev_dict['u'] = out_pu
        prev_dict['v'] = out_pv

def lamb_vortex_merge_python(vormat, mincoreratio, maxcoreratio, radiusratio, maxcoresizeratio):
    """For CCSVM merging
    
    Suggested parameters 0.9, 1.1, 0.5, 2.0
    
    """
    mergecount=0
    cntbfor = len(vormat['x'])
    # try merge! WARNING NOT EQUIPED TO DEAL WITH VORTICES OF
    # NEGATIVE ROTATION!
    imerge=0
    total_vorts = len(vormat['x'])
    cndodgy = 0

    # faster to not call item everytime in loop!!!
    vormat_x = vormat['x']
    vormat_y = vormat['y']
    vormat_s = vormat['s']
    vormat_c = vormat['c']

    while imerge < total_vorts:
        #print imerge, '/', len(vormat_x)
        mergelist = []
        merge_x = []
        merge_y = []
        merge_s = []
        merge_c = []
        #cumulativestrength=0;
        # still potential left to find merges begining search
        thiscore = vormat_c[imerge]
        mincore = mincoreratio * thiscore
        maxcore = maxcoreratio * thiscore
        maxcoresize = thiscore * maxcoresizeratio
        FAKEFUSERADIUS = thiscore * radiusratio
        FAKEFUSERADIUSsqrd = FAKEFUSERADIUS ** 2
        # searching for available merges
        for jmerge in xrange(imerge + 1, total_vorts):
            # only add to list if core size is within min<o<max ratios...(suggest 0.9 + 1.1)
            if mincore <= vormat_c[jmerge] <= maxcore: 
            # setting merge radius at some distance this current particle would split! maybe should set it on external field factors...ie merge at a radius such that it will not split next step from diffusion/travel!!.
                rsqrd = (vormat_x[imerge] - vormat_x[jmerge]) ** 2 + (vormat_y[imerge] - vormat_y[jmerge]) ** 2
                if rsqrd < FAKEFUSERADIUSsqrd:
                    # !!!!!!USING MARKS HODGE-PODGE MERGING
                    # APPROXIMATIONS...THESE ARE possibly NOT VALID TO ROSSI ALGORITHM.
                    #if rsqrd <= R*minsigmasqrd
                    #cumulativestrength = cumulativestrength + vormat(jmerge,3);
                    #if cumulativestrength < mergingepsilon
                    mergelist.append(jmerge)
                    merge_x.append(vormat_x[jmerge])
                    merge_y.append(vormat_y[jmerge])
                    merge_s.append(vormat_s[jmerge])
                    merge_c.append(vormat_c[jmerge])
                    mergecount += 1

        if len(mergelist) > 0:
            # Merges to this particle are possible so doing them now
            #DONT FORGET TO INCLUDE ITSELF IN THE MERGE!!!
            mergelist.insert(0, imerge)
            merge_x.insert(0, vormat_x[imerge])
            merge_y.insert(0, vormat_y[imerge])
            merge_s.insert(0, vormat_s[imerge])
            merge_c.insert(0, vormat_c[imerge])

            mergecount += 1
            # Get the strength of the new vortex, conserving 0th moment
            stn = sum(merge_s)
            # Get the x-position of the new vortex conserving 1st moment
            xn = sum((merge_s[i] * merge_x[i] for i in xrange(len(mergelist)))) / stn
            # Get the y-position of the new <F12>vortex conserving 1st
            # moment
            yn = sum((merge_s[i] * merge_y[i] for i in xrange(len(mergelist)))) / stn
            # Get the radial distances to the old cores
            Rs = [(yn - merge_y[i]) ** 2 + (xn - merge_x[i]) ** 2 for i in xrange(len(mergelist))]
            # Get the core-size of the new vortex
            # cn = sqrt(sum(strg.*(4*(cors.^2) + Rs))/(4*stn));         # This is from Rossi
            cn = math.sqrt(sum((merge_s[i] * (merge_c[i]**2 + Rs[i]) for i in xrange(len(mergelist))))/stn)          # This is modified for my core (no 4)  !!!!!!! ABSOLUTE VALUE!!!!!!

            # Impose a limit on the maximum core size that can be created
            if cn > maxcoresize:                                                       
            #!!!!!!!! CAREFUL HERE TOO !!!!!!!!!!!!
                cn = maxcoresize
                cndodgy += 1

            # modifying vortex in current position
            vormat_x[imerge] = xn
            vormat_y[imerge] = yn
            vormat_s[imerge] = stn
            vormat_c[imerge] = cn
            
            mergelist.reverse()
            # only works in reverse as indices would be rooted otherwise
            for i in mergelist: 
                # deleting merged vortices except for the original which should be at the back
                if i != imerge:
                    del vormat_x[i]
                    del vormat_y[i]
                    del vormat_s[i]
                    del vormat_c[i]
            
            total_vorts = len(vormat_x)

        # updating for next loop
        imerge=imerge + 1
    
    hassplit = 0
    
    if cndodgy > 0:
        print('doing ' + str(cndodgy) + ' cn = maxcoresize dodgies')

    print 'merged', mergecount,'vortices.', cntbfor-len(vormat_x),'vortices deleted'
    

    ## FINISHED MERGING VORTICES IF NECESSARY ---------------------
    
    
def real_main():
    #test_lamb_vortex_split()
    #test_lamb_vortex_merge()
    #test_lamb_vortex_merge(True)
    #test_lamb_vortex_combination()
    test_lamb_vortex_prev_u()

def test_lamb_vortex_split():
    # test alpha = 0,0.5,1
    # create vormats
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat_keys = ['x', 'y', 's', 'c']
    split_inds = [0]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 1
    num_children = 4
    print '\n','num_children', num_children, ', split_alpha', split_alpha, 'and vormat before :::'
    print [vormat[0], vormat[1],vormat[2],vormat[3]]
    print '==========--------========AFTER:----->'
    newvormat = lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print [newvormat[0], newvormat[1],newvormat[2],newvormat[3]]
    
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 0
    num_children = 4
    print '\n','\n','num_children', num_children, ', split_alpha', split_alpha, 'and vormat before :::'
    print [vormat[0], vormat[1],vormat[2],vormat[3]]
    print '==========--------========AFTER:----->'
    newvormat = lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print [newvormat[0], newvormat[1],newvormat[2],newvormat[3]]
    
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 0.5
    num_children = 4
    print '\n','num_children', num_children, ', split_alpha', split_alpha, 'and vormat before :::'
    print [vormat[0], vormat[1],vormat[2],vormat[3]]
    print '==========--------========AFTER:----->'
    newvormat = lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print [newvormat[0], newvormat[1],newvormat[2],newvormat[3]]
    
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 0.5
    num_children = 6
    print '\n','num_children', num_children, ', split_alpha', split_alpha, 'and vormat before :::'
    print [vormat[0], vormat[1],vormat[2],vormat[3]]
    print '==========--------========AFTER:----->'
    newvormat = lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print [newvormat[0], newvormat[1],newvormat[2],newvormat[3]]
    print "validated successfully 2/12/2009"


def test_lamb_vortex_merge(fmm=False):
    vormat_keys = ['x', 'y', 's', 'c']
    print '\n\n \n \n'
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.9
    maxo = 1.1
    cutrad = 1.0
    maxcutrad = 2
    print '\n \n \t test single element'
    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print vormat,
    print '==========--------========AFTER:----->'
    if not fmm:
        newvormat = lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        newvormat = lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)

    print vormat,

    
    x = [0.0]*3
    y = [0.0]*3
    s = [1.0]*3
    c = [1.0]*3
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    vormat['s'][0] *= -1
    mino = 0.9
    maxo = 1.1
    cutrad = 1.0
    maxcutrad = 2
    print '\n \n \t two elements on top another + one negative'
    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print vormat,
    print '==========--------========AFTER:----->'
    if not fmm:
        lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        newvormat = lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)
    print vormat,


    x = [-2, -3, -1, -2]
    y = [0.0] * 3 + [1.11]
    s = [1.0] * 4
    c = [1.1] * 4
    
    x += [2, 1, 3, 1] 
    y += [0.0] * 3 + [1.11]
    s += [1.0] * 4
    c += [1.1] * 4
    
    x += [0.0] 
    y += [0.0]
    s += [1.0]
    c += [1.1]

    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.9
    maxo = 1.1
    cutrad = 1.0
    maxcutrad = 2
    
    print '\n \n \t 4 elements, seperated by axis, and offset, + 2 a little far'
    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print vormat,
    print '==========--------========AFTER:----->'
    if not fmm:
        newvormat = lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        newvormat = lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)
    print vormat,

    print '\n \n \t testing conservation of str and 1st moment and 2nd moment. in a random bunch of elements...maybe do a split if need be'
    numtotest = 500
    x = [random.random() for i in xrange(numtotest)]
    y = [random.random() for i in xrange(numtotest)]
    s = [random.random() for i in xrange(numtotest)]
    c = [random.random() for i in xrange(numtotest)]

    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.9
    maxo = 1.1
    cutrad = 0.5
    maxcutrad = 2

    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print '1st, 2nd moment respectively;', vormat.moment(0), vormat.moment(1)
    print '==========--------========AFTER:----->'
    if not fmm:
        newvormat = lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        newvormat = lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)
    print '1st, 2nd moment respectively;', vormat.moment(0), vormat.moment(1)
    print "validated successfully 9/12/2009"

def test_lamb_vortex_prev_u():
    # create vormats
    x = [0.0]
    y = [0.0]
    s = [1.0]
    c = [1.0]
    vmt = [x, y, s, c]
    vormat_keys = ['x', 'y', 's', 'c']
    split_inds = [0]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 0.85
    num_children = 4
    print '\n','no prev_u',
    print '==========--------========AFTER:----->'
    lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print [vormat[0], vormat[1], vormat[2], vormat[3]]
    
    x = [0.0,1.0]
    y = [0.0,0.0]
    s = [1.0,1.0]
    c = [1.0,1.0]
    prev_u = [33.0, 33.7]
    prev_v = [12.0, 33.7]
    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    split_alpha = 0.85
    num_children = 4
    split_inds = [0,1]
    print '\n','\n','with prev_u and prev_v of before :::'
    print [vormat[0], vormat[1],vormat[2],vormat[3], prev_u, prev_v]
    print '==========--------========AFTER:----->'
    lamb_vortex_split(vormat, split_inds, num_children, split_alpha, prev_vormat_u=prev_u, prev_vormat_v=prev_v)
    newvormat = vormat
    print [newvormat[0], newvormat[1],newvormat[2],newvormat[3], prev_u, prev_v]
    
    print "split prev_u validated successfully 15/01/2010"

    x = [-2, -3, -1, -2]
    y = [0.0] * 3 + [1.11]
    s = [1.0] * 4
    c = [1.1] * 4
    
    x += [2, 1, 3, 1] 
    y += [0.0] * 3 + [1.11]
    s += [1.0] * 4
    c += [1.1] * 4
    
    x += [0.0] 
    y += [0.0]
    s += [1.0]
    c += [1.1]

    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.9
    maxo = 1.1
    cutrad = 1.0
    maxcutrad = 2
    fmm=True
    print '!!!!!TESTING NO PREV_u'
    print '\n \n \t 4 elements, seperated by axis, and offset, + 2 a little far'
    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print vormat,
    print '==========--------========AFTER:----->'
    if not fmm:
        lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)
    print vormat,

    
    
    x = [-2, -3, -1, -2]
    y = [0.0] * 3 + [1.11]
    s = [1.0] * 4
    c = [1.1] * 4
    
    x += [2, 1, 3, 1] 
    y += [0.0] * 3 + [1.11]
    s += [1.0] * 4
    c += [1.1] * 4
    
    x += [0.0] 
    y += [0.0]
    s += [1.0]
    c += [1.1]

    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.9
    maxo = 1.1
    cutrad = 1.0
    maxcutrad = 2
    fmm=True
    pu = range(9)
    pv = range(9)
    pu[0] = -1
    pv[0] = -1
    prev_dict = {'u':pu, 'v':pv}
    print '\n \n \t 4 elements, seperated by axis, and offset, + 2 a little far'
    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print vormat,
    print pu
    print pv
    print '==========--------========AFTER:----->'
    if not fmm:
        lamb_vortex_merge_python(vormat, mino,maxo,cutrad, maxcutrad)
    else:
        lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad, prev_dict=prev_dict)
        pu = prev_dict['u']
        pv = prev_dict['v']
    print vormat,
    print pu
    print pv

    print '\n \n \t comparing prev_u and prev_v'
    print "validated successfully 15/01/2010"

def test_lamb_vortex_combination():
    vormat_keys = ['x', 'y', 's', 'c']
    print '\n \n \t test combination, random particles, split + merge, conservation of 0,1,2 moments'
    numtotest = 20000
    x = [random.random() for i in xrange(numtotest)]
    y = [random.random() for i in xrange(numtotest)]
    s = [random.random() for i in xrange(numtotest)]
    c = [random.random() for i in xrange(numtotest)]

    vmt = [x, y, s, c]
    vormat = Elemat(vmt, vormat_keys)
    mino = 0.85
    maxo = 1.15
    cutrad = 0.25
    maxcutrad = 2
    split_alpha = 0.75
    num_children = 4
    split_inds = range(numtotest)

    print 'radius of cuttof rswtio', cutrad, ', and vormat before :::'
    print '1st, 2nd moment respectively;', vormat.moment(0), vormat.moment(1)
    print '==========--------========AFTER SPLIT:----->'
    newvormat = lamb_vortex_split(vormat, split_inds, num_children, split_alpha)
    print '1st, 2nd moment respectively;', vormat.moment(0), vormat.moment(1)
    print '==========--------========AFTER MERGE:----->'
    newvormat = lamb_vortex_merge(vormat, mino,maxo,cutrad, maxcutrad)
    print '1st, 2nd moment respectively;', vormat.moment(0), vormat.moment(1)
    
    print "validated successfully 9/12/2009"


def profile_main():
    """profiler function"""

    import cProfile, pstats
    prof = cProfile.Profile()
    prof = prof.runctx("real_main()", globals(), locals())
    stats = pstats.Stats(prof)
    stats.sort_stats("cumulative")  # Or cumulative
    stats.print_stats(80)  # 80 = how many to print
    # The rest is optional.
    # stats.print_callees()
    # stats.print_callers()

if __name__ == '__main__':
    real_main()
