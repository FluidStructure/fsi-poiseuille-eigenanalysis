from numpy import mat, reshape, linspace, reshape, append, random, ones, delete
from numpy import meshgrid, bmat, array, size, zeros, transpose, cumsum, sign, mod
from scipy.io import savemat, loadmat
import numpy as np
import os
import pdb

class finiteDifference():

    def getInputShape(self,n):
        if size(n) == size(n,0):
            nr = size(n,0)
            nc = 0
            R = zeros(nr)
        else:
            nr = size(n,0)
            nc = size(n,1)
            R = zeros((nr,nc))
        return nr, R    

    def d4dx(self,n,N,dx):
        # Calculates the 4th order derivative of vector 'n' using finite difference
        # method applied to every 'N'th element with an element spacing of 'dx'
        # NOTE: if 'v' is a square Identity matrix of type 'array' then the output is a square matrix
        # that includes the finite difference coefficients
        (nn,R) = self.getInputShape(n)
        # Only hinged-hinged conditions for now
        i = 0;R[i] = 0.0 - 4.0*0.0 + 0.0*n[i] - 0.0*n[i+N] + 0.0*n[i+2*N]
        i = 1;R[i] = 0.0 - 4.0*n[i-N] + 5.0*n[i] - 4.0*n[i+N] + n[i+2*N]
        i = nn-2;R[i] = n[i-2*N] - 4.0*n[i-N] + 5.0*n[i] - 4.0*n[i+N] + 0.0
        i = nn-1;R[i] = 0.0*n[i-2*N] - 0.0*n[i-N] + 0.0*n[i] - 0.0 + 0.0
        for i in range(2,nn-2):
            R[i] = n[i-2*N] - 4.0*n[i-N] + 6.0*n[i] - 4.0*n[i+N] + n[i+2*N]
        R = R/(dx**4.0)
        return R

    def d2dx(self,v,N,dx):
        (nn,R) = self.getInputShape(v)
        dxsq = dx**2.0
        for i in range(0,nn):
            if (i < N):
                R[i] = (0 - 2.0*v[i] + 1.0*v[i+N])      # Hinged end condition (upstream virtual node = 0)
            elif (i >= nn-N):
                R[i] = (2.0*v[i-N] - 2.0*v[i] + 0.0)       # Free end conditions (downstream gradient = 0)
            else:
                R[i] = (v[i-N] - 2.0*v[i] + v[i+N])
        R = R/(dx**2.0)
        return R

    def d1dx(self,v,N,dx,order=2):
        # Evaluate the x-direction gradient
        (nn,R) = self.getInputShape(v)
        if order == 1:
            # Use a first order (upwind) scheme
            for i in range(0,nn):
                if (i < N):
                    R[i] = v[i]/dx        # Hinged end condition (virtual node = 0)
                else:
                    R[i] = (v[i] - v[i-N])/dx
        elif order == 2:
            # Use a second order upwinding scheme
            for i in range(0,nn):
                if i< N:
                    #R[i] = (-1*0 +4*0 - 3*v[i])/(2*dx)
                    R[i] = v[i]/dx              # First order upwind on the first node
                elif i < 2*N:
                    #R[i] = (-1*0 +4*v[i-N] - 3*v[i])/(2*dx)
                    R[i] = (v[i] - v[i-N])/dx   # First order upwind on the second node
                else:
                    R[i] = (1.0*v[i-2*N] - 4.0*v[i-N] + 3.0*v[i])/(2.0*dx)
        elif order == 3:
            # Use a third order upwinding scheme
            for i in range(0,nn):
                if i< N:
                    # First order for the first node
                    #R[i] = (-1*0 + 6*0 - 3*v[i] -2*v[i+N])/(6*dx)
                    R[i] = v[i]/dx
                elif i < 2*N:
                    R[i] = (1.0*0 - 6.0*v[i-N] + 3.0*v[i] + 2.0*v[i+N])/(6.0*dx)
                elif i >= nn-N:
                    # Second order upwind for the second last node
                    R[i] = (1.0*v[i-2*N] - 4.0*v[i-N] + 3.0*v[i])/(2.0*dx)
                else:
                    R[i] = (1.0*v[i-2*N] - 6.0*v[i-N] + 3.0*v[i] + 2.0*v[i+N])/(6.0*dx)
        return R


class chebychev():

    def chebdif(self,N,M):

        from numpy import eye, floor, ceil, sin, cos, ones
        from numpy import zeros, matrix, mat, linspace, transpose
        from numpy import pi, arange, repeat, dot, flipud, fliplr
        from numpy import bmat, power, diag, sum, logical_not
        from numpy import multiply, array, size
        from scipy.linalg.basic import toeplitz

        """
        function [x, DM] = chebdif(N, M)

        %  The function [x, DM] =  chebdif(N,M) computes the differentiation 
        %  matrices D1, D2, ..., DM on Chebyshev nodes. 
        % 
        %  Input:
        %  N:        Size of differentiation matrix.        
        %  M:        Number of derivatives required (integer).
        %  Note:     0 < M <= N-1.
        %
        %  Output:
        %  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
        %
        %  The code implements two strategies for enhanced 
        %  accuracy suggested by W. Don and S. Solomonoff in 
        %  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
        %  The two strategies are (a) the use of trigonometric 
        %  identities to avoid the computation of differences 
        %  x(k)-x(j) and (b) the use of the "flipping trick"
        %  which is necessary since sin t can be computed to high
        %  relative precision when t is small whereas sin (pi-t) cannot.
        %  Note added May 2003:  It may, in fact, be slightly better not to
        %  implement the strategies (a) and (b).   Please consult the following
        %  paper for details:   "Spectral Differencing with a Twist", by
        %  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 

        %  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by 
        %  JACW, May 2003.
        """
        
        I = eye(N)                          # Identity matrix.     
        L = logical_not(1 - I)              #L = logical(I) # Logical identity matrix.

        n1 = floor(N/2) 
        n2  = ceil(N/2)                     # Indices used for flipping trick.

        k = matrix(linspace(0,N-1,N))
        k = transpose(k)                    # k = [0:N-1]' # Compute theta vector.
        
        th = k*pi/(N-1)

        x = sin(pi*transpose(mat(arange(N-1,1-N-1,-2)))/(2*(N-1)))    # Compute Chebyshev points.

        T = repeat(th/2,N,axis=1)              
        DX = multiply(2*sin(transpose(T)+T),sin(transpose(T)-T))        # Trigonometric identity. 
        DX = bmat([[DX[0:n1,:]],[-1*flipud(fliplr(DX[0:n2,:]))]])       # Flipping trick. 
        DX[L] = ones(N)                                                 # Put 1's on the main diagonal of DX.
        
        C = mat(toeplitz(power(-1,k)))               # C is the matrix with 
        C[0,:] = C[0,:]*2
        C[N-1,:] = C[N-1,:]*2                   # entries c(k)/c(j)
        C[:,0] = C[:,0]/2
        C[:,N-1] = C[:,N-1]/2
        
        Z = 1/DX                                # Z contains entries 1/(x(k)-x(j))
        Z[L] = zeros(N)                         # with zeros on the diagonal.
        
        D = eye(N)                              # D contains diff. matrices.
        
        DM = {}                                 # Initialize a dictionary for holding D matrices 
        
        for ell in range(1,M+1):
            D = multiply((ell*Z),multiply(C,repeat(transpose(mat(diag(D))),N,axis=1))-D)      # Off-diagonals
            D[L] = -1*sum(transpose(array(D)),axis=0)                                     # Correct main diagonal of D
            DM[ell] = D                                                     # Store current D in DM

        return x, DM


    def cheb2bc(self,N,g):
    
        from numpy import eye, floor, ceil, sin, cos, ones
        from numpy import zeros, matrix, mat, linspace, transpose
        from numpy import pi, arange, repeat, dot, flipud, fliplr
        from numpy import bmat, power, diag, sum, logical_not
        from numpy import multiply, array, size
        from scipy.linalg.basic import toeplitz    
    
        """
        function [xt,D2t,D1t,phip,phim]=cheb2bc(N,g);
        %
        % Program for computing first and second derivative matrices and
        % and boundary condition functions for 2 point boundary conditions
        %
        %  a_1 u(1)  + b_1 u'(1)  = c_1
        %  a_N u(-1) + b_N u'(-1) = c_N
        %
        %
        % INPUT 
        % N        =  number of Chebyshev points in [-1,1]
        % g        =  boundary condition matrix = [a_1 b_1 c_1; a_N b_N c_N]
        % 
        % OUTPUT  
        % xt       =  Chebyshev points corresponding to rows and columns
        %             of D1t and D2t
        % D1t      =  1st derivative matrix incorporating bc
        % D2t      =  2nd derivative matrix incorporating bc
        % phip     =  1st and 2nd derivative of bc function at x=1
        %             (array with 2 columns)
        % phim     =  1st and 2nd derivative of bc function at x=-1 
        %             (array with 2 columns)

        % S.C. Reddy, J.A.C. Weideman  1998
        """


        # Get differentiation matrices
        (x,DM) = self.chebdif(N,2)
        D0=mat(eye(N,N))
        D1=DM[1]
        D2=DM[2]

        # extract boundary condition coefficients
        a1=g[0,0]; b1=g[0,1]; c1=g[0,2];
        aN=g[1,0]; bN=g[1,1]; cN=g[1,2];

        # Case 0: Invalid boundary condition information
        if ((a1==0.0 and b1==0.0) or (aN==0.0 and bN==0.0)):
            print 'Invalid boundary condition information (no output)'
            
        elif ((b1==0.0) and (bN==0.0)):                 # Dirichlet/Dirichlet 
            
            J=mat(arange(2,N)) - 1
            K=transpose(mat(arange(2,N))) - 1
            D1t=D1[J,K]
            D2t=D2[J,K]
            phip=c1*bmat([D1[K,0], D2[K,0]])/a1     # phi_+
            phim=cN*bmat([D1[K,N], D2[K,N]])/aN     # phi_- 
            xt=x[K]                                 # node vector 
            
        elif ((b1!=0.0) and (bN==0.0)):             # Dirichlet x=-1, Robin x=1
            
            J = mat(arange(2,N)) - 1
            K = transpose(mat(arange(1,N))) - 1
            xjrow=2*power(sin(J*pi/2/(N-1)),2)      # 1-x_j, using trig identity
            xkcol=2*power(sin(K*pi/2/(N-1)),2)      # 1-x_k, using trig identity
            oner=mat(ones((size(xkcol,axis=0),size(xkcol,axis=1))))     # column of ones

            fac0 = oner*(1/xjrow)                   # matrix -1/(1-x_j)
            fac1 = xkcol*(1/xjrow)                  # matrix (1-x_k)/(1-x_j)
            D1t = multiply(fac1,D1[K,J])-multiply(fac0,D0[K,J])
            D2t = multiply(fac1,D2[K,J])-2*multiply(fac0,D1[K,J]) 

            cfac = D1[0,0]+a1/b1                    # compute phi'_1, phi''_1
            fcol1 = -cfac*D0[K,0]+multiply((1+cfac*xkcol),D1[K,0])
            fcol2 = -2*cfac*D1[K,0]+multiply((1+cfac*xkcol),D2[K,0])
            D1t  = bmat([fcol1,D1t])
            D2t  = bmat([fcol2,D2t])

            phim = multiply(xkcol,(D1(K,N-1)/2))-(D0(K,N-1)/2)      # phi'_-, phi''_- 
            phim = cN*bmat([phim,multiply(xkcol,D2(K,N-1)/2)-D1(K,N-1)])/aN

            phip = multiply(-1*xkcol,D1(K,0))+D0(K,0)           # phi'_+, phi''_+ 
            phip = c1*bmat([phip,-1*multiply(xkcol,D2(K,0))+2*D1(K,0)])/b1
            
            xt = x[K]                                           # node vector

        elif ((b1==0.0) and (bN!=0.0)):
            # Case 3: Dirichlet at x=1 and Neumann or Robin boundary x=-1.

            J = mat(arange(2,N)) - 1 
            K = transpose(mat(arange(1,N+1))) - 1
            xjrow=2*power(cos(J*pi/2/(N-1)),2)      # 1+x_j, using trig identity
            xkcol=2*power(cos(K*pi/2/(N-1)),2)      # 1+x_k, using trig identity
            oner=mat(ones((size(xkcol,axis=0),size(xkcol,axis=1))))                # column of ones

            fac0 = oner*(1/xjrow)                   # matrix 1/(1+x_j)
            fac1 = xkcol*(1/xjrow)                  # matrix (1+x_k)/(1+x_j)
            D1t = multiply(fac1,D1[K,J])+multiply(fac0,D0[K,J])
            D2t = multiply(fac1,D2[K,J])+2*multiply(fac0,D1[K,J]) 

            cfac = D1[N-1,N-1]+aN/bN                # compute phi'_N, phi''_N
            lcol1 = -cfac*D0[K,N-1]+multiply((1-cfac*xkcol),D1[K,N-1])
            lcol2 = -2*cfac*D1[K,N-1]+multiply((1-cfac*xkcol),D2[K,N-1])
            D1t  = bmat([D1t,lcol1])
            D2t  = bmat([D2t,lcol2])

            phim = multiply(xkcol,(D1[K,N-1]))+(D0[K,N-1])      # phi'_-, phi''_- 
            phim = cN*bmat([phim,multiply(xkcol,D2[K,N-1])+2*D1[K,N-1]])/bN

            phip = multiply(xkcol,(D1[K,0]/2))+D0[K,0]          # phi'_+, phi''_+ 
            phip = c1*bmat([phip,multiply(xkcol,D2[K,0]/2)+D1[K,0]])/a1

            xt = x[K];                                          # node vector

        elif ((b1!=0.0) and (bN!=0.0)):

            # Case 4: Neumann or Robin boundary conditions at both endpoints. 

            J = mat(arange(2,N)) - 1
            K = transpose(mat(arange(1,N+1))) - 1
            xkcol0=power(sin((K)*pi/(N-1)),2)           # 1-x_k^2 using trig identity
            xkcol1=transpose(-2*x[K])                   # -2*x_k 
            xkcol2=-2*mat(ones((size(xkcol0,axis=0),size(xkcol0,axis=1))))  # -2
            xjrow=1/power(sin((J)*pi/(N-1)),2)          # 1-x_j^2 using trig identity

            fac0=xkcol0*xjrow
            fac1=xkcol1*xjrow
            fac2=xkcol2*xjrow

            D1t=multiply(fac0,D1[K,J])+multiply(fac1,D0[K,J])
            D2t=multiply(fac0,D2[K,J])+multiply(2*fac1,D1[K,J])+multiply(fac2,D0[K,J])

            omx=power(sin((K)*pi/2/(N-1)),2)                # (1-x_k)/2 
            opx=power(cos((K)*pi/2/(N-1)),2)                # (1+x_k)/2

            r0=opx+(0.5+D1[0,0]+a1/b1)*xkcol0/2             # compute phi'_1, phi''_1
            r1=0.5-(0.5+D1[0,0]+a1/b1)*x
            r2=-0.5-D1[0,0]-a1/b1
            rcol1=multiply(r0,D1[K,0])+multiply(r1,D0[K,0])
            rcol2=multiply(r0,D2[K,0])+multiply(2*r1,D1[K,0])+multiply(r2,D0[K,0])

            l0=omx+(0.5-D1[N-1,N-1]-aN/bN)*xkcol0/2            # compute phi'_N, phi''_N
            l1=-0.5+(D1[N-1,N-1]+aN/bN-0.5)*x
            l2=D1[N-1,N-1]+aN/bN-0.5
            lcol1=multiply(l0,D1[K,N-1])+multiply(l1,D0[K,N-1])
            lcol2=multiply(l0,D2[K,N-1])+2*multiply(l1,D1[K,N-1])+multiply(l2,D0[K,N-1])

            D1t=bmat([rcol1,D1t,lcol1])
            D2t=bmat([rcol2,D2t,lcol2])

            phim1=(multiply(xkcol0,D1[K,N-1])+multiply(xkcol1,D0[K,N-1]))/2
            phim2=(multiply(xkcol0,D2[K,N-1])+2*multiply(xkcol1,D1[K,N-1])+multiply(xkcol2,D0[K,N-1]))/2
            phim=cN*bmat([phim1,phim2])/bN                 # compute phi'_-, phi''_-

            phip1=(multiply(-xkcol0,D1[K,0])-multiply(xkcol1,D0[K,0]))/2
            phip2=(multiply(-xkcol0,D2[K,0])-2*multiply(xkcol1,D1[K,0])-multiply(xkcol2,D0[K,0]))/2
            phip=c1*bmat([phip1,phip2])/b1                 # compute phi'_+, phi''_+

            xt=x[K]                                 # node vector
            
        return xt,D2t,D1t,phip,phim

class ode():

    def __init__(self):
        self.cc = 0
        self.plotHalf = True
        self.outPath = ''

    def demoOde45(self):
        #
        # Test of the ode45 function using a predator-prey model
        #

        def tfun(t,y):
            print("t=%f"%t)
            yd = [0.0]*2
            yd[0] =-.5*y[0] + .01*y[0]*y[1];
            yd[1] = .5*y[1] -.01*y[0]*y[1];
            return yd

        tslot=[0.0,50.0];
        #tslot=[i/4.0 for i in range(201)]
        yinit = [80.0,100.0];

        #o = ode()
        #o.ode45(tfun,tslot,yinit)
        self.ode45(tfun,tslot,yinit)

        ui = raw_input("Press ENTER to end...")

    def odeSaveMat(self,t,y):
        print 'OdeSaveMat, t = ' + str(t) + ', cc = ' + str(self.cc)
        outDict = {}
        outDict['t'] = t
        outDict['y'] = y
        # Check if path exists and make it if not
        if not os.path.exists(self.outPath):
            os.makedirs(self.outPath)
        fn = 'TStep_%04d.mat'%self.cc
        fn = self.outPath + fn
        savemat(fn,outDict)
        self.cc += 1
        return None
        
        
    def odeplot(self,t,y):
        #print 'OdePlot, t = ' + str(t) + ', y = ' + str(y)
        print 'OdePlot, t = ' + str(t) + ', cc = ' + str(self.cc)
        if self.plotHalf == True:
            y = [y[i] for i in xrange(len(y)/2)]
        if self.cc == 0:
            fig1 = plt.figure(figsize=(16,4))
            plt.ion()
            plt.grid(True)
            #plt.hold(True)
            if len(y) > 4:
                plt.plot(y)
                self.ylims = plt.ylim()
            else:
                self.t = [t]
                self.y = [];self.y.append(y)
                plt.plot(self.t,self.y)
        else:
            plt.clf()
            if len(y) > 4:
                plt.plot(y)
                plt.ylim(self.ylims)
            else:
                self.t.append(t)
                self.y.append(y)
                plt.plot(self.t,self.y)
        plt.draw()
        self.cc += 1
        #ui = raw_input("Press ENTER to end...")
        return None

    def ode45(self, vfun, vslot, vinit, RelTol=None, AbsTol=None,
        NormControl=True, MaxStep=None, NonNegative=None, Mass=None,
        InitialStep=None, Refine=0, OutputFcn=None,
        OutputSel=None, OutputSave=-1, vparams=None):
        """
        %# Copyright (C) 2006-2008, Thomas Treichl <treichl@users.sourceforge.net>
        %# OdePkg - A package for solving ordinary differential equations and more
        %#
        %# This program is free software; you can redistribute it and/or modify
        %# it under the terms of the GNU General Public License as published by
        %# the Free Software Foundation; either version 2 of the License, or
        %# (at your option) any later version.
        %#
        %# This program is distributed in the hope that it will be useful,
        %# but WITHOUT ANY WARRANTY; without even the implied warranty of
        %# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        %# GNU General Public License for more details.
        %#
        %# You should have received a copy of the GNU General Public License
        %# along with this program; If not, see <http://www.gnu.org/licenses/>.

        %# -*- texinfo -*-
        %# @deftypefn  {Function File} {[@var{}] =} ode45 (@var{@@fun}, @var{slot}, @var{init}, [@var{opt}], [@var{par1}, @var{par2}, @dots{}])
        %# @deftypefnx {Command} {[@var{sol}] =} ode45 (@var{@@fun}, @var{slot}, @var{init}, [@var{opt}], [@var{par1}, @var{par2}, @dots{}])
        %# @deftypefnx {Command} {[@var{t}, @var{y}, [@var{xe}, @var{ye}, @var{ie}]] =} ode45 (@var{@@fun}, @var{slot}, @var{init}, [@var{opt}], [@var{par1}, @var{par2}, @dots{}])
        %#
        %# This function file can be used to solve a set of non--stiff ordinary differential equations (non--stiff ODEs) or non--stiff differential algebraic equations (non--stiff DAEs) with the well known explicit Runge--Kutta method of order (4,5).
        %#
        %# If this function is called with no return argument then plot the solution over time in a figure window while solving the set of ODEs that are defined in a function and specified by the function handle @var{@@fun}. The second input argument @var{slot} is a double vector that defines the time slot, @var{init} is a double vector that defines the initial values of the states, @var{opt} can optionally be a structure array that keeps the options created with the command @command{odeset} and @var{par1}, @var{par2}, @dots{} can optionally be other input arguments of any type that have to be passed to the function defined by @var{@@fun}.
        %#
        %# If this function is called with one return argument then return the solution @var{sol} of type structure array after solving the set of ODEs. The solution @var{sol} has the fields @var{x} of type double column vector for the steps chosen by the solver, @var{y} of type double column vector for the solutions at each time step of @var{x}, @var{solver} of type string for the solver name and optionally the extended time stamp information @var{xe}, the extended solution information @var{ye} and the extended index information @var{ie} all of type double column vector that keep the informations of the event function if an event function handle is set in the option argument @var{opt}.
        %#
        %# If this function is called with more than one return argument then return the time stamps @var{t}, the solution values @var{y} and optionally the extended time stamp information @var{xe}, the extended solution information @var{ye} and the extended index information @var{ie} all of type double column vector.
        %#
        %# For example, solve an anonymous implementation of the Van der Pol equation
        %#
        %# @example
        %# fvdb = @@(vt,vy) [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
        %#
        %# vopt = odeset ("RelTol", 1e-3, "AbsTol", 1e-3, \
        %#          "NormControl", "on", "OutputFcn", @@odeplot);
        %# ode45 (fvdb, [0 20], [2 0], vopt);
        %# @end example
        %# @end deftypefn
        %#
        %# @seealso{odepkg}

        %# ChangeLog:
        %#   20010703 the function file "ode45.m" was written by Marc Compere
        %#     under the GPL for the use with this software. This function has been
        %#     taken as a base for the following implementation.
        %#   20060810, Thomas Treichl
        %#     This function was adapted to the new syntax that is used by the
        %#     new OdePkg for Octave and is compatible to Matlab's ode45.
        """
        
        eps = np.finfo(float).eps
        
        if (len(vslot) > 2):     # Step size checking
            vstepsizefixed = True;
        else:
            vstepsizefixed = False;
        
        if (RelTol == None and not vstepsizefixed):
            RelTol = 1e-6;
            print('OdePkg:InvalidArgument, Option "RelTol" not set, new value %f is used'%RelTol);
        elif (RelTol != None and vstepsizefixed):
            print('OdePkg:InvalidArgument, Option "RelTol" will be ignored if fixed time stamps are given');
        
        if (AbsTol == None and not vstepsizefixed):
            AbsTol = 1e-6
            print('OdePkg:InvalidArgument, Option "AbsTol" not set, new value %f is used'%AbsTol);
        elif (AbsTol != None and vstepsizefixed):
            print('OdePkg:InvalidArgument, Option "AbsTol" will be ignored if fixed time stamps are given');

        if NormControl == True:
            vnormcontrol = True
        else:
            vnormcontrol = False
        
        if not NonNegative == None:
            if Mass == None:
                vhavenonnegative = True;
            else:
                vhavenonnegative = False;
                print('OdePkg:InvalidArgument, Option "NonNegative" will be ignored if mass matrix is set');
        else:
            vhavenonnegative = False

        if (OutputFcn==None):
            OutputFcn = self.odeSaveMat;
            vhaveoutputfunction = True;
        else:
            vhaveoutputfunction = True;

        if OutputSel==None:
            vhaveoutputselection = False;
        else:
            vhaveoutputselection = True;

        if (Refine > 0.0):
            vhaverefine = True
        else:
            vhaverefine = False
        
        if (InitialStep == None) and not vstepsizefixed:
            InitialStep = (vslot[1] - vslot[0]) / 10.0;
            InitialStep = InitialStep / (10.0 ** Refine);
            print('OdePkg:InvalidArgument, Option "InitialStep" not set, new value %f is used'%InitialStep);

        if (MaxStep == None) and not vstepsizefixed:
            MaxStep = (vslot[1] - vslot[0]) / 10.0;
            print('OdePkg:InvalidArgument, Option "MaxStep" not set, new value %f is used'%MaxStep);

        if (Mass != None) and (type(Mass)==matrix):
            vhavemasshandle = False; vmass = Mass;      # constant mass
        elif (Mass != None):
            vhavemasshandle = True;                     # mass defined by a function handle
        else: # no mass matrix - creating a diag-matrix of ones for mass
            vhavemasshandle = False;                    # vmass = diag (ones (length (vinit), 1), 0);

        if (OutputSel == None):
            vhaveoutputselection = False;
        else:
            vhaveoutputselection = True;
        
        
        # Starting the initialisation of the core solver ode45 
        vtimestamp  = vslot[0];             # timestamp = start time
        vtimelength = len(vslot);           # length needed if fixed steps
        vtimestop   = vslot[vtimelength-1]; # stop time = last value
        vdirection  = sign (vtimestop);     # Flag for direction to solve

        if not vstepsizefixed:
            vstepsize = InitialStep;
            vminstepsize = (vtimestop - vtimestamp) / (1.0/eps);
        else: # If step size is given then use the fixed time steps
            vstepsize = vslot[1] - vslot[0];
            vminstepsize = sign (vstepsize) * eps;
        
        vretvaltime = []
        vretvalresult = []
        
        vretvaltime.append(vtimestamp);   # first timestamp output
        if (vhaveoutputselection):  # first solution output
            vretvalresult.append(vinit(OutputSel))
        else:
            vretvalresult.append(vinit)

        vpow = 1.0/5.0;                                         # 20071016, reported by Luis Randez
        va = mat([[0.0, 0.0, 0.0, 0.0, 0.0],                    # The Runge-Kutta-Fehlberg 4(5) coefficients
            [1.0/4.0, 0.0, 0.0, 0.0, 0.0],                      # Coefficients proved on 20060827
            [3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0],                # See p.91 in Ascher & Petzold
            [1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0],
            [439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0],
            [-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0]])
        # 4th and 5th order b-coefficients
        vb4 = mat([25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0]).transpose();
        vb5 = mat([16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0]).transpose();
        vc = [sum(va[i]) for i in xrange(len(va))];


        # The solver main loop - stop if the endpoint has been reached
        vcntloop = 2; vcntcycles = 1; vu = mat(vinit); 
        vk = vu.transpose() * zeros((1,6));
        #vk = mat([[0.0]*6 for i in xrange(len(vu))]);
        vcntiter = 0; vunhandledtermination = True; vcntsave = 2;
        while ((vdirection * (vtimestamp) < vdirection * (vtimestop)) and (vdirection * (vstepsize) >= vdirection * (vminstepsize))):
            # Hit the endpoint of the time slot exactely
            if ((vtimestamp + vstepsize) > vdirection * vtimestop):
                vstepsize = vtimestop - vdirection * vtimestamp;

            # Estimate the six results when using this solver
            for j in xrange(6):
                vthetime  = vtimestamp + vc[j] * vstepsize;
                vtheinput = vu.transpose() + vstepsize * vk[:,0:j] * va[j,0:j].transpose();
                if (vhavemasshandle):                       # Handle only the dynamic mass matrix,
                    if (vmassdependence):                   # constant mass matrices have already
                        vmass = Mass(vthetime,vtheinput)    # been set before (if any)
                    else:                                   # if (vmassdependence == false)
                        vmass = Mass(vthetime)              # then we only have the time argument
                    vk[:,j] = np.linalg.solve(vmass,vfun(vthetime, vtheinput));
                else:
                    tin = [vtheinput[i,0] for i in xrange(size(vtheinput))]
                    vk[:,j] = mat(vfun(vthetime, tin)).transpose();

            # Compute the 4th and the 5th order estimation
            y4 = vu.transpose() + vstepsize * (vk * vb4);
            y5 = vu.transpose() + vstepsize * (vk * vb5);
            if (vhavenonnegative):
              vu[NonNegative] = abs (vu[NonNegative]);
              y4[NonNegative] = abs (y4[NonNegative]);
              y5[NonNegative] = abs (y5[NonNegative]);
            vSaveVUForRefine = vu;

            # Calculate the absolute local truncation error and the acceptable error
            if not vstepsizefixed:
                if not vnormcontrol:
                    vdelta = array(abs (y5 - y4));
                    vtau = max(max(RelTol * abs (vu.transpose())), AbsTol);
                else:
                    vdelta = array([np.linalg.norm (y5 - y4, ord=np.inf)]);
                    vtau = max (RelTol * max (np.linalg.norm (vu.transpose(), ord=np.inf), 1.0), AbsTol);
            else:   # if (vstepsizefixed == true)
                vdelta = array([1.0]); vtau = 2;

            # If the error is acceptable then update the vretval variables
            if (np.all (vdelta <= vtau)):
                vtimestamp = vtimestamp + vstepsize;
                vu = y5.transpose();          # MC2001: the higher order estimation as "local extrapolation"
                # Save the solution every vodeoptions.OutputSave steps             
                if (mod (vcntloop-1,OutputSave) == 0) and (OutputSave >= 0):
                    if (vhaveoutputselection):
                        vretvaltime.append(vtimestamp);
                        vretvalresult.append(vu[OutputSel]);
                    else:
                        vretvaltime.append(vtimestamp);
                        vretvalresult.append(vu);
                    vcntsave = vcntsave + 1;    
                vcntloop = vcntloop + 1; vcntiter = 0;

                # Call plot only if a valid result has been found, therefore this
                # code fragment has moved here. Stop integration if plot function
                # returns false
                if (vhaveoutputfunction):
                    for vcnt in xrange(Refine+1):           # Approximation between told and t
                        if (vhaverefine):                   # Do interpolation
                            vapproxtime = (vcnt + 1) * vstepsize / (Refine + 2);
                            vapproxvals = vSaveVUForRefine.transpose() + vapproxtime * (vk * vb5);
                            vapproxtime = (vtimestamp - vstepsize) + vapproxtime;
                        else:
                            vapproxvals = vu.transpose();
                            vapproxtime = vtimestamp;
                        if (vhaveoutputselection):
                            vapproxvals = vapproxvals(OutputSel);
                        tin = [vapproxvals[i,0] for i in xrange(size(vapproxvals))]
                        vpltret = OutputFcn(vapproxtime, tin);          
                        if vpltret:     # Leave refinement loop
                            break
    ##                    if (vpltret):   # Leave main loop
    ##                        vunhandledtermination = False;
    ##                        break;

            # Update the step size for the next integration step
            if not vstepsizefixed:
                # 20080425, reported by Marco Caliari
                # vdelta cannot be negative (because of the absolute value that
                # has been introduced) but it could be 0, then replace the zeros 
                # with the maximum value of vdelta
                #print(vdelta)
                for i in xrange(len(vdelta)):
                    if vdelta[i] == 0.0:
                        vdelta[i] = max(vdelta)
                # It could happen that max (vdelta) == 0 (ie. that the original
                # vdelta was 0), in that case we double the previous vstepsize
                for i in xrange(len(vdelta)):
                    if vdelta[i] == 0.0:
                        vdelta[i] = multiply(vtau,(0.4 ** (1.0 / vpow)));

                if (vdirection == 1):
                    vstepsize = min (MaxStep, min (0.8 * vstepsize * (vtau / vdelta) ** vpow));
                else:
                    vstepsize = max (MaxStep, max (0.8 * vstepsize * (vtau / vdelta) ** vpow));

            else:   # if (vstepsizefixed)
                if (vcntloop <= vtimelength):
                    vstepsize = vslot[vcntloop-1] - vslot[vcntloop-2];
                else: # Get out of the main integration loop
                    break;

            # Update counters that count the number of iteration cycles
            vcntcycles = vcntcycles + 1;    # Needed for cost statistics
            vcntiter = vcntiter + 1;        # Needed to find iteration problems

            # Stop solving because the last 1000 steps no successful valid
            # value has been found
            if (vcntiter >= 5000):
              merr("Solving has not been successful")