#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Library for calculating theoretical boundary layer profiles.

AUTHOR - Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
MODIFIED - 2009.08

"""

def calc_lam_bound_lay_vel(
    bl_y,
    bnd_lay_height=1,
    umean=1,
    mode=None,
    ):
    """#function [outVel,vort] = calcLamBoundLayVel(blY,bndLayHeight,Umean)
    
    function [outVel,vort] = calcLamBoundLayVel(blY,bndLayHeight,Umean)
    This is a function that returns the velocity profile
    of a Blasius boundary layer over a flat plate using
    the Pohlhausen 4th order approximation to achieve
    the solution.
    
    FORMAT:
    [udivU,vort] = pohl(n)
    INPUTS:
    n -         is the normalised y-coordinate with respect to the
    boundary layer thickness
    OUTPUTS:
    udivU -    is the normalised velocity profile with respect
    to the mean flow velocity
    udivU = u/U
    vort -    is the normalised vorticity with respect to the
    boundary layer thickness and Umean
    normvort = -1*vort*blthick/(2*Umean)
    
    ##outVel = Umean * sqrt(blY./bndLayHeight);
    
    ## Set points that are above \"n=1\" and less than
    ## \"n=0\" equal to \"0\"
    
    """

    if hasattr(bl_y, '__iter__'):

    # input and output an array

        out_vel = []
        vort = []
        for height in bl_y:
            (svel, svort) = calc_lam_bound_lay_vel_single(height,
                    bnd_lay_height, umean)
            out_vel.append(svel)
            vort.append(svort)
    else:
        (out_vel, vort) = calc_lam_bound_lay_vel_single(bl_y,
                bnd_lay_height, umean)

    if not mode is None:
        out_vel = vort
    return out_vel


def calc_lam_bound_lay_vel_single(height, bnd_lay_height, umean):
    """#function [outVel,vort] = calcLamBoundLayVel(blY,bndLayHeight,Umean)
    
    function [outVel,vort] = calcLamBoundLayVel(blY,bndLayHeight,Umean)
    This is a function that returns the velocity profile
    of a Blasius boundary layer over a flat plate using
    the Pohlhausen 4th order approximation to achieve
    the solution.
    
    FORMAT:
    [udivU,vort] = pohl(n)
    INPUTS:
    n -         is the normalised y-coordinate with respect to the
    boundary layer thickness
    OUTPUTS:
    udivU -    is the normalised velocity profile with respect
    to the mean flow velocity
    udivU = u/U
    vort -    is the normalised vorticity with respect to the
    boundary layer thickness and Umean
    normvort = -1*vort*blthick/(2*Umean)
    
    ##outVel = Umean * sqrt(blY./bndLayHeight);
    
    ## Set points that are above \"n=1\" and less than
    ## \"n=0\" equal to \"0\"
    
    """

    normalized_height = height / bnd_lay_height
    udiv_u = 0.0
    if normalized_height > 1.0:
        udiv_u = 1.0
    elif normalized_height >= 0.0:
        # Calculate the  velocity profile
        udiv_u = 2.0 * normalized_height - 2.0 * normalized_height ** 3\
             + normalized_height ** 4
    else: 
        raise ValueError('normalized_height < 0.0...no such thing')

    out_vel = udiv_u * umean

    # Calculate the vorticity profile achieved through simple differentiation

    vort = ((2 * umean) * (1 - 3 * normalized_height ** 2 + 2
             * normalized_height ** 3)) / bnd_lay_height
    return (out_vel, vort)


if __name__ == '__main__':

    # print calcLamBoundLayVel.__doc__

    print '######### IN TESTING'
    C = [(k - 1) * 0.1 for k in range(13)]
    (A, B) = calc_lam_bound_lay_vel(C, 1, 1)
    print A
    print '             '
    print B
