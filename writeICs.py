#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This code uses Jarrads FMM to construct the large square influence coefficient
matrices for fluid-fluid, fluid-wall, wall-fluid and wall-wall for later use
in some linear-algebra method such as "ode45" or "eigs"
"""

from numpy import mat, reshape, linspace, reshape, append, random, ones
from numpy import meshgrid, bmat, array, size, zeros, transpose
from scipy.io import savemat, loadmat
import velocitylib
import pdb

def fluidICs(xcF,ycF,xcW,ycW,dx):
    
    # Get the number of fluid elementsfrom numpy import mat, reshape, linspace, reshape, append, random, ones
    Nf = len(xcF)
    Nw = len(xcW)
    
    # Generate the list of all the elements (fluid - top wall - bottom wall)
    X = xcF + xcW
    Y = ycF + ycW
    
    # Solve the velocity at all elements due to fluid elements
    INff = zeros((Nf,Nf))
    ITff = zeros((Nf,Nf))
    INfw = zeros((Nw,Nf))
    ITfw = zeros((Nw,Nf))
    for i in xrange(Nf):
        print 'Calculating fluid-all influence for fluid element ' + str(i) + ' of ' + str(Nf)
        x = xcF[i]
        y = ycF[i]
        (U,V) = velocitylib.vortex_element_vel([x-(dx/2.0)],[y],[1.0],[x+(dx/2.0)],[y],[1.0],X,Y)
        INff[:,i] = V[0:Nf]
        ITff[:,i] = U[0:Nf]
        INfw[:,i] = V[Nf:Nf+(2*Nw)]
        ITfw[:,i] = U[Nf:Nf+(2*Nw)]

    return INff,ITff,INfw,ITfw
    
    
def wallICs(xcF,ycF,xcW,ycW,dx):
    
    # Get the number of fluid elementsfrom numpy import mat, reshape, linspace, reshape, append, random, ones
    Nf = len(xcF)
    Nw = len(xcW)
    
    # Generate the list of all the elements (fluid - top wall - bottom wall)
    X = xcF + xcW
    Y = ycF + ycW
    
    # Solve the velocity at all elements due to wall elements
    INww = zeros((Nw,Nw))
    ITww = zeros((Nw,Nw))
    INwf = zeros((Nf,Nw))
    ITwf = zeros((Nf,Nw))
    for i in xrange(Nw):
        print 'Calculating wall-all influence for wall element ' + str(i) + ' of ' + str(Nw)
        x = xcW[i]
        y = ycW[i]
        (U,V) = velocitylib.source_element_vel([x-(dx/2.0)],[y],[1.0],[x+(dx/2.0)],[y],[1.0],X,Y)
        INwf[:,i] = V[0:Nf]
        ITwf[:,i] = U[0:Nf]
        INww[:,i] = V[Nf:Nf+(2*Nw)]
        ITww[:,i] = U[Nf:Nf+(2*Nw)]
        # Correct the self-inflence (diagonals) for wall-wall influence coefficients
        if y > 0:
            INww[i,i] = -0.5

    return INww,ITww,INwf,ITwf

if __name__ == "__main__":
    print 'Constructing square influence-coefficent matrices...'

    # Load in variables from .mat file
    invec = loadmat('VARS.mat',struct_as_record=True)
    # Extract important variables
    xcF = [invec["pcxxv"][i][0] for i in xrange(len(invec["pcxxv"]))]
    ycF = [invec["pcyyv"][i][0] for i in xrange(len(invec["pcyyv"]))]
    dx = invec["dxc"][0][0]
    xc = [invec["pcx"][0][i] for i in xrange(len(invec["pcx"][0]))]
    xcW = xc + xc
    ycW = [1.0]*len(xc) + [-1.0]*len(xc)
    
    (INff,ITff,INfw,ITfw) = fluidICs(xcF,ycF,xcW,ycW,dx)
    (INww,ITww,INwf,ITwf) = wallICs(xcF,ycF,xcW,ycW,dx)
    
    # Save the variables to a .mat file (matlab/octave format)
    outvec = {}
    outvec["INff"] = INff
    outvec["ITff"] = ITff
    outvec["INfw"] = INfw
    outvec["ITfw"] = ITfw
    outvec["INww"] = INww
    outvec["ITww"] = ITww
    outvec["INwf"] = INwf
    outvec["ITwf"] = ITwf
    print 'Writing influence coefficient matrices to ICs.mat'
    savemat('ICs.mat',outvec)
    print 'Finished...'