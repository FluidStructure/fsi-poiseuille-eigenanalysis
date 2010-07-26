from scipy import reshape, real, imag, array, exp, pi, linspace, divide, multiply, iscomplexobj
from scipy import append, size, concatenate, dot, eye, mod, linalg, sum, column_stack
from scipy.io import savemat, loadmat
from scipy.sparse.linalg import gmres, minres, bicgstab, LinearOperator
from scipy.sparse.linalg.eigen.arpack.speigs import ARPACK_eigs
from scipy.sparse.linalg.eigen.arpack import eigen
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib.pyplot import imshow, axis, contour, contourf
from numpy import mat, reshape, linspace, reshape, append, random, ones, delete
from numpy import meshgrid, bmat, array, size, zeros, transpose, cumsum, add
import numpy as np
import os, sys, subprocess

import pdb
import cProfile, pstats

import sys
sys.path.append('../python_tools.git/')
import velocitylib
import octFuncs

def merr(self,msg):
    raise TypeError(msg)

class parameters():

    def parseBoolString(self,s):
        if s == 'True' or s == 'true':
            b = True
        elif s == 'False' or s == 'false':
            b = False
        return b

    def parseParameterFile(self,FILE='PARAMETERS.dat'):

        # Read the parameters for the system
        f = open(FILE, 'r')
        text = f.read()
        f.close()

        lines = text.splitlines()

        p = []
        for line in lines:
            value = line.split('#')[0]
            value = value.replace(' ','')
            if not value == '':
                p.append(value)
        
        return p
    
    class simulation():

        def __init__(self):
            c = parameters()
            self.parseBoolString = c.parseBoolString
            p = c.parseParameterFile(FILE='PARAMETERS.dat')
            self.setParameters(p)

        def setParameters(self,p):

            self.R = float(p[0])
            self.chebN = int(p[1])
            self.Nup = int(p[2])
            self.Nco = int(p[3])
            self.Ndn = int(p[4])
            self.LT = float(p[5])
            self.rhow = float(p[6])
            self.hw = float(p[7])
            self.d = float(p[8])
            self.E = float(p[9])
            self.K = float(p[10])
            self.fluidOnly = self.parseBoolString(p[11])
            self.allowSlip = self.parseBoolString(p[12])
            self.deterministicBCs = self.parseBoolString(p[13])
            self.periodicBCs = self.parseBoolString(p[14])
            self.Ndomains = int(p[15])
            self.method = p[16]
            self.Nthreads = int(p[17])
            self.octave = p[18]
            self.solver = p[19]
            self.Nsteps = float(p[20])
            self.maxStep = float(p[21])
            self.invTol = float(p[22])
            self.eigTol = float(p[23])
            self.apl = float(p[24])
            
            # Some calculated fixed variables for the simulation
            self.Nx = self.Nup+self.Nco+self.Ndn                # Total number of panels along the length of the domain
            self.Nf = self.chebN*self.Nx                        # Total number of fluid elements
            self.dx = self.LT/self.Nx                           # Element length in the x-direction
            self.B = self.E*(self.hw**3.0)/12.0                 # Flexural Rigidity of the wall
            self.FS = self.LT*(self.Nup/self.Nx)                # Starting coordinate of the flexible wall section
            self.L = self.LT*(float(self.Nco)/float(self.Nx))   # Length of the wall compliant (m)
            delta = 1.0                                         # Convert Reynolds number to a boundary layer thickness for numerical calculation (delta = half channel width)
            ymax = 2.0;                                         # For truncation of the computational domain in the wall-normal y-domain.  The maximum y-value
            Uinf = 1.0;                                         # Mean flow velocity (m/s)
            rho = 1.0;                                          # Flow density (m/s)
            self.mu = Uinf*rho*delta/self.R;                    # Flow viscosity
            self.nu = self.mu/rho;                              # Flow kinematic viscosity
    
    class postProcessing():

        def __init__(self):
            c = parameters()
            self.parseBoolString = c.parseBoolString
            p = c.parseParameterFile(FILE='POSTPRPARS.dat')
            self.setParameters(p)
        
        class general():
            pass
        class movie():
            pass
        def setParameters(self,p):
            g = self.general
            g.modeNum = int(p[0])                                   # Mode number to plot
            g.resultsDir = p[1]                                     # Path to results directory
            g.normalizeSpatialGrowth = self.parseBoolString(p[2])   # Normalises spatial growth (to see absolute growth/decay)
            g.ignoreTemporalGrowth = self.parseBoolString(p[3])     # Useful for visualisation
            g.maxFactor = float(p[4])                               # Plot contours from +/-maxFactor*max(abs(V))

            m = self.movie
            m.P = int(p[5])                                         # Number of time periods to plot movie over
            m.Np = int(p[6])                                        # Number of time steps per period
            m.includeEndTime = self.parseBoolString(p[7])           # Do not include endtime if movie is going to loop
            m.fps = int(p[8])                                       # Number of frames per second
            m.vcodec = p[9]                                         # Video codec to use

class geometry():
    
    def __init__(self):
        p = parameters.simulation()
        self.makeGeometry(p)
    
    def makeGeometry(self,p):

        c = octFuncs.chebychev()

        # Get the chebyshev grid points and 1st and 2nd order differentiation matrices over [-1,1]
        bcsy = mat([[0,1,0],[0,1,0]])           # Boundary conditions
        (ynbc,d2dy,d1dy,phip,phim) = c.cheb2bc(p.chebN,bcsy);
        self.d2dy = d2dy
        
        # Get the node points of each cell and cell widths
        yn = [0.0]*(p.chebN+1)
        yn[0] = float(ynbc[0])
        for ss in range(0,p.chebN-1):
            yn[ss+1] = (float(ynbc[ss+1]) + float(ynbc[ss]))/2.0
        yn[p.chebN] = float(ynbc[p.chebN-1])
        self.yn = yn
        
        # Get width and centres of each cell
        self.dy = [0.0]*(p.chebN)
        for s in range(0,p.chebN):
            self.dy[s] = abs(yn[s+1] - yn[s])
        
        # Get the position of the flow elements (collocation points)
        pcy = ynbc
        pcy[0] = (yn[0] + yn[1])/2
        pcy[p.chebN-1] = (yn[p.chebN] + yn[p.chebN-1])/2.0
        self.yc = [float(pcy[i]) for i in xrange(len(pcy))]
        #pcy = transpose(pcy[1:size(pcy)-1])
        
        # Get the x-positions of the flow elements
        self.xn = linspace(0.0,p.LT,p.Nx+1)
        pcx = linspace(0.5*p.dx,p.LT-(0.5*p.dx),p.Nx)
        
        # Get the x-coordinates of the nodes of the cmopliant wall only
        FS = p.LT*(p.Nup/p.Nx); 	    # Starting point of the flexible wall section
        L = p.LT*(p.Nco/p.Nx);    	    # Length of the wall compliant (m)
        self.xnCW = linspace(FS,FS+L,(p.Nco+1))
        
        # Create a meshgrid and then reshape to get vector of fluid element positions
        (xcF,ycF) = meshgrid(pcx,pcy)
        self.xcF = reshape(xcF,(p.chebN)*p.Nx,1)
        self.ycF = reshape(ycF,(p.chebN)*p.Nx,1)
        
        # Create a vector of wall element positions
        self.xc = pcx
        self.xcW = append(pcx,pcx)
        self.ycW = append(((pcx*0.0)+1.0),((pcx*0.0)-1.0))

class fmmMethod():


    def __init__(self):
        self.parameters = parameters.simulation()
        self.geometry = geometry()
        #------------------
        # Get some critical geometry that is used often
        g = self.geometry
        p = self.parameters
        self.xo = None
        self.eCount = 0
        self.eVec = []
        self.x0 = None
        # Adjust panel tolerance on FMM call (and warn if necessary)
        self.panTol = 0.001
        if ((1.0-g.yc[0])/p.dx) < 0.001:
            self.panTol = ((1.0-g.yc[0])/p.dx)/2.0
            print('ADJUSTING PANEL TOLERANCE FOR FMM TO ' + str(self.panTol))
        # Get near-wall fluid element geometry
        (g.xcFnw,g.xcFfw) = self.getNearAndFarWallElements(g.xcF)
        (g.ycFnw,g.ycFfw) = self.getNearAndFarWallElements(g.ycF)
        # Generate periodic elements
        (g.XLfw,g.XRfw) = self.genPeriodicElements(g.xcFfw)
        (g.XLnw,g.XRnw) = self.genPeriodicElements(g.xcFnw)
        (g.XLw,g.XRw) = self.genPeriodicElements(g.xcW)

    def getNearAndFarWallElements(self,x):
        p = self.parameters
        xnw = [x[i*p.chebN] for i in xrange(p.Nx)] + [x[((i+1)*(p.chebN))-1] for i in xrange(p.Nx)]
        xfw = []
        for i in xrange(p.Nx):
            xfw += list(x[(i*p.chebN + 1):((i+1)*p.chebN - 1)])
        return xnw, xfw

    def runEigSolver(self):
        p = self.parameters
        g = self.geometry

        Ny = p.chebN
        Nx = p.Nx
        if p.deterministicBCs == True:
            Ny = Ny - 2

        if p.fluidOnly == False:
            N = Ny*Nx + (p.Nco-1)*2
        else:
            N = Ny*Nx
        
        def tfun(y):
            self.eCount += 1    # Update the iteration counter
            self.eVec = y       # Update the initial guess for the future
            print('Calling EIGEN tfun.  Iteration = ' + str(self.eCount))
            if iscomplexobj(y):
                print("---WARNING!  Input vector is complex!")
                merr('Input vector is complex!!')
                #yd_r = self.multRHS(real(y))
                #yd_i = self.multRHS(imag(y))
                #yd = [complex(yd_r[i],yd_i[i]) for i in xrange(len(y))]
            else:
                if p.fluidOnly == True:
                    yd = self.multRHS(y)
                else:
                    RHSvector = self.multRHS(y)
                    LHS = LinearOperator( (N,N), matvec=self.multLHS, dtype='float64' )
                    (yd,F) = bicgstab(LHS,transpose(mat(RHSvector)),M=pCond,tol=p.invTol,x0=RHSvector)
                    if F != 0:
                        merr('Iterative matrix inverse did not converge.')
                    self.x0 = yd    # Update initial guess
            return yd

        ## Using the eigen call
        RHS = LinearOperator( (N,N), matvec=tfun, dtype='float64' )
        #print(int(p.Nsteps))
        #(e,v) = eigen(RHS, k=6, M=None, sigma=None, which='LR', v0=None, ncv=None, maxiter=int(p.Nsteps), tol=p.eigTol, return_eigenvectors=True)
        (e,v) = eigen(RHS, k=6, M=None, sigma=None, which='LR', v0=None, ncv=None, maxiter=None, tol=p.eigTol, return_eigenvectors=True)
        ## Using ARPACK_eigs
        #(e,v) = ARPACK_eigs(tfun, N, 6, which='LR', ncv=None, tol=p.eigTol)
        print(e)
        
        # Make the output filename
        fname = self.makeEigFilename()
        fname = 'results/' + fname
        
        # Save the result
        outDict = {}
        outDict['evals'] = e
        outDict['Veigs'] = v
        outDict['y0'] = self.eVec
        savemat(fname,outDict)

    def makeEigFilename(self):
        p = self.parameters
        g = self.geometry
        fname = 'eigs_'
        if p.fluidOnly == True:
            fname += 'FLO_'
        else:
            fname += 'FSI_'
        fname = fname + str(p.chebN) + 'x' + str(p.Nx)
        fname = fname + '_R' + str(int(p.R))
        fname = fname + '.mat'
        return fname

    def runOdeSolver(self):
        p = self.parameters
        g = self.geometry

        Ny = p.chebN
        Nx = p.Nx
        if p.deterministicBCs == True:
            Ny = Ny - 2

        if p.fluidOnly == False:
            N = Ny*Nx + (p.Nco-1)*2
        else:
            N = Ny*Nx

        def tfun(t,y):
            if p.fluidOnly == True:
                yd = self.multRHS(y)
            else:
                RHSvector = self.multRHS(y)
                LHS = LinearOperator( (N,N), matvec=self.multLHS, dtype='float64' )
                (yd,F) = bicgstab(LHS,transpose(mat(RHSvector)),tol=p.invTol,x0=RHSvector)
                if F != 0:
                    merr('Iterative matrix inverse did not converge.')
            return yd
        yinit = [0.0]*N
##        # Initialize with a spot disturbance
##        nn = int(round(Nx/4.0)*Ny + round(Ny/2.0))
##        yinit[nn] = 1e-2;yinit[nn+Ny*20] = -1e-2
##        yinit[nn-1] = 1e-2;yinit[nn-1+Ny*20] = -1e-2
        # Initialize with a sine-wave disturbance at the centreline
        NyHalf = int(round(Ny/2.0))
        for i in xrange(Nx):
            yinit[i*Ny + NyHalf] = 1e-2*np.sin(2.0*pi*g.xc[i]/(p.LT/10.0))
            yinit[i*Ny + NyHalf-1] = 1e-2*np.sin(2.0*pi*g.xc[i]/(p.LT/10.0))
        # Set the timeslots
        tslot = [0,p.Nsteps]
        # Call the ode45 solver
        o = octFuncs.ode()
        path = 'results/' + str(int(p.R))
        if p.fluidOnly == True:
            path = path + '/fluidOnly/'
        else:
            path = path + '/FSI/'
        o.outPath = path
        o.ode45(tfun,tslot,yinit,MaxStep=p.maxStep,InitialStep=p.maxStep)

    def callOctave(self):
        # Save generalised matrix to A.mat
        outDict = {}
        outDict['A'] = self.RHSGeneral
        savemat('A.mat',outDict)
        # Run octave script from command line
        os.system(self.parameters.octave + ' --eval callOctaveSolvers.m')

    def getFlowAndWallParts(self,v,Nf,Nn):
        vf = v[0:Nf]
        if Nf+2*Nn == len(v):
            vw = v[Nf:Nf+2*Nn]
        else:
            vw = []
        return vf,vw

    def getElementNumbers(self):
        p = self.parameters
        g = self.geometry
        
        Nn = p.Nco - 1          # The number of nodes in the compliant panel section
        if p.deterministicBCs == True:
            Ny = p.chebN - 2        # The number of vertical fluid elements
        else:
            merr('SORRY!  Can only use deterministic boundary conditions in FMM method for now')
        
        # Get fluid and wall components of input vector "v"
        Nx = p.Nx
        Nf = Ny*p.Nx
        return Nx, Ny, Nf, Nn

    def multRHS(self,v):
        print('-Multiplying by RHS matrix')
        p = self.parameters
        g = self.geometry
        f = octFuncs.finiteDifference()
        # Ensure that the input vector is a list
        v = [vv for vv in v]
        # The the number of various elements and flow and wall parts
        (Nx, Ny, Nf, Nn) = self.getElementNumbers()
        (vf,vw) = self.getFlowAndWallParts(v,Nf,Nn)
        
        # Multiply flow element strengths "v" by element height to get singularity strengths
        dy = g.dy[1:p.chebN-1]
        vfs = [dy[mod(i,Ny)]*vf[i] for i in xrange(Nf)]
        
        # Solve the velocity at all elements due to fluid elements (excluding nearest-to-wall fluid elements)
        Np = len(g.XLfw)/len(g.xcFfw)
        (Ufa,Vfa) = velocitylib.vortex_element_vel(g.XLfw,g.ycFfw*Np,vfs*Np,g.XRfw,g.ycFfw*Np,vfs*Np,(g.xcFfw+list(g.xcW)),(g.ycFfw+list(g.ycW)),threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
        Uf = Ufa[:Nf];Ufw = Ufa[Nf:]
        Vf = Vfa[:Nf];Vfw = Vfa[Nf:]
       
        # Get the total normal flux through wall elements
        Vw = Vfw
        if p.fluidOnly == False:
            Vnodes = [0.0] + vw[0:p.Nco-1] + [0.0]
            for c in range(p.Nco):
                i = p.Nx+p.Nup+c
                Vw[i] -= (Vnodes[c] + Vnodes[c+1])/2.0

        # Solve for wall source strengths and nearest-to-wall flow element strengths together
        RHSw = append(multiply(Vw,-1.0),multiply(Ufw,-1.0))
        sigma = self.getWallElementStrengths(RHSw)

        # Get the normal velocity at fluid elements due to wall source and near-wall vortex elements and
        # Add to normal velocity at fluid elements due to themselves (already calculated)
        vNWs = [s*g.dy[0] for s in sigma[(Nx*2):(Nx*4)]]
        vNW = [s for s in sigma[(Nx*2):(Nx*4)]]
        (Uwf,Vwf) = velocitylib.source_element_vel(g.XLw,list(g.ycW)*Np,sigma[0:2*Nx]*Np,g.XRw,list(g.ycW)*Np,sigma[0:2*Nx]*Np,g.xcFfw,g.ycFfw,threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
        Vf = [Vf[i] + Vwf[i] for i in xrange(len(Vwf))]
        (Uwf,Vwf) = velocitylib.vortex_element_vel(g.XLnw,g.ycFnw*Np,vNWs*Np,g.XRnw,g.ycFnw*Np,vNWs*Np,g.xcFfw,g.ycFfw,threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
        Vf = [Vf[i] + Vwf[i] for i in xrange(len(Vwf))]
        
        # Calculate convective components (finite difference) for RHS of fluid transport
        if p.periodicBCs == 'False':
            # Just do a normal finite diference
            d1dx = f.d1dx(vf,Ny,p.dx,order=2)
            d2dx = f.d2dx(vf,Ny,p.dx)
        else:
            # Do a cyclic finite difference
            vfc = vf[-2*Ny:] + vf + vf[:2*Ny]
            #-----------
            d1dx = list(f.d1dx(vfc,Ny,p.dx,order=2))
            d2dx = list(f.d2dx(vfc,Ny,p.dx))
            del d1dx[-2*Ny:]
            del d2dx[-2*Ny:]
            del d1dx[:2*Ny]
            del d2dx[:2*Ny]
        # Get chebyshev second-order gradients in the y-direction (must add in near-wall vortex strengths)
        vfIncNearWall = []
        #pdb.set_trace()
        ##vNW = [0.0]*len(g.xcFnw)    # For debugging only
        for i in xrange(p.Nx):
            vfIncNearWall += [vNW[i]]
            vfIncNearWall += vf[i*Ny:(i+1)*Ny]
            vfIncNearWall += [vNW[i+p.Nx]]
        d2dy = self.d2dyCheb(vfIncNearWall,g.d2dy);
        (d2dynw,d2dy) = self.getNearAndFarWallElements(list(d2dy))
        
        # Put together the RIGHT HAND SIDE of the FLUID PART
        vdot = [0.0]*Nf
        for i in xrange(Nf):
            # Get local mean flow
            U = -1.0*(1 - (g.ycFfw[i]**2.0))
            # Get vdot
            #vdot[i] = U*d1dx[i]
            #vdot[i] = U*d1dx[i] + p.nu*(d2dx[i] + d2dy[i])
            vdot[i] = Vf[i]*2.0 + U*d1dx[i] + p.nu*(d2dx[i] + d2dy[i])
        
        # Add the wall terms
        if p.fluidOnly == False:
            vdotWall = [0.0]*len(vw)
            v1 = vw[:Nn]
            v2 = vw[Nn:]
            # v1_dot = v2
            for i in xrange(Nn):
                vdotWall[i] += v2[i]
            # v2_dot = RHS of wall equation
            d4dx = f.d4dx([0.0]+v1+[0.0],1,p.dx)
            d4dx = d4dx[1:len(d4dx)-1]
            for i in xrange(Nn):
                vdotWall[i+Nn] += -1.0*(p.B*d4dx[i] + p.K)*v1[i] - p.d*v2[i]
            # Append this onto the output vector
            vdot += vdotWall
        
        return vdot
        
    def multLHS(self,v):
        print('-Multiplying by LHS matrix')
        p = self.parameters
        g = self.geometry
        f = octFuncs.finiteDifference()
        # Ensure that the input vector is a list
        v = [vv for vv in v]
        # The the number of various elements and flow and wall parts
        (Nx, Ny, Nf, Nn) = self.getElementNumbers()
        (vf,vw) = self.getFlowAndWallParts(v,Nf,Nn)
        v1 = vw[:Nn]
        v2 = vw[Nn:]

        # Multiply flow element strengths "v" by element height to get singularity strengths
        dy = g.dy[1:p.chebN-1]
        vfs = [dy[mod(i,Ny)]*vf[i] for i in xrange(Nf)]

        # Ones for all the fluid elements
        vdot = [vv for vv in vf]

        # Ones for v1 = v2 Part of the wall equations
        vdot += [vv for vv in v1]
        
        # v2 = RHS of wall equation

        # Solve the velocity at all elements due to fluid elements (excluding nearest-to-wall fluid elements)
        Np = len(g.XLfw)/len(g.xcFfw)
        (Ufw,Vfw) = velocitylib.vortex_element_vel(g.XLfw,g.ycFfw*Np,vfs*Np,g.XRfw,g.ycFfw*Np,vfs*Np,list(g.xcW),list(g.ycW),threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
       
        # Get the total normal flux through wall elements
        Vw = Vfw
        if p.fluidOnly == False:
            for c in range(p.Nco):
                i = p.Nx+p.Nup+c
                Vw[i] -= (vw[c] + vw[c+1])/2.0

        # Solve for wall source strengths and nearest-to-wall flow element strengths together
        RHSw = append(multiply(Vw,-1.0),multiply(Ufw,-1.0))
        sigma = self.getWallElementStrengths(RHSw)
       
        sVLW = sigma[Nx*3:Nx*4]                         # Get the strength of the vertices along the lower wall
        PF = [0.0] + [g.dy[0]*p.dx*s for s in sVLW]
        PF = cumsum(PF)
        PF = [-1.0*pp for pp in PF]
        Pav = sum(PF)/len(PF)
        PFp = [(pp-Pav) for pp in PF]
        
##        PFp = [0.0] + [-1.0*g.dy[0]*p.dx*(sVLW[i]+sVLW[i+1])/2.0 for i in xrange(len(sVLW)-1)] + [0.0]
        
        # Add pressures and wall density etc to the LHS
        vdot += [v2[i]*p.rhow*p.hw + PFp[p.Nup+1+i] for i in xrange(len(v2))]

        return vdot

    def getWallElementStrengths(self,RHSv):
        p = self.parameters
        g = self.geometry
        (Nx, Ny, Nf, Nn) = self.getElementNumbers()
        # Solve for wall source strengths and nearest-to-wall flow element strengths together
        def preCond(xx):
            xo = xx
            xo[0*Nx:1*Nx] = (xo[0*Nx:1*Nx])*-2.0                  # Upper source elements
            xo[1*Nx:2*Nx] = (xo[1*Nx:2*Nx])*2.0                   # Lower source elements
            xo[2*Nx:3*Nx] = (xo[2*Nx:3*Nx])*2.0*(1/g.dy[0])       # Upper vortex elements
            xo[3*Nx:4*Nx] = (xo[3*Nx:4*Nx])*-2.0*(1/g.dy[0])      # Lower vortex elements
            return xo
        def multINTww(xx):
            #print("Doing INTww*x...")
            Np = len(g.XLnw)/len(g.ycFnw)   # Number of periodic domains (for copying y-coordinates)
            
            # Evaulate the matrix-vector product of source/sink influence coefficients using Jarrads FMM
            #print('Evaluating matrix-vector product: INTww*X...')
            ss = [x for x in xx[0:(Nx*2)]]
            (us,vs) = velocitylib.source_element_vel(g.XLw,list(g.ycW)*Np,ss*Np,g.XRw,list(g.ycW)*Np,ss*Np,list(g.xcW),list(g.ycW),threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
            # Adjust the upper panels to have a self-influence of -0.5
            for i in range(0,Nx):
                vs[i] = vs[i] - ss[i]
                
            # Evaluate the matrix-vector product of wall-vortex influence coefficients
            sv = [x*g.dy[0] for x in xx[(Nx*2):(Nx*4)]]
            (uv,vv) = velocitylib.vortex_element_vel(g.XLnw,g.ycFnw*Np,sv*Np,g.XRnw,g.ycFnw*Np,sv*Np,list(g.xcW),list(g.ycW),threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
            # Construct the output vector
            ov = [vs[i] + vv[i] for i in xrange(len(vv))] + [us[i] + uv[i] for i in xrange(len(uv))]
            return ov
        self.xo = preCond(RHSv)
        INTww = LinearOperator( (Nx*2*2,Nx*2*2), matvec=multINTww, dtype='float64' )
        pCond = LinearOperator( (Nx*2*2,Nx*2*2), matvec=preCond, dtype='float64' )
        #(sigma,F) = minres(INTww,transpose(mat(RHSv)),M=pCond,tol=p.invTol,x0=self.xo)
        (sigma,F) = gmres(INTww,transpose(mat(RHSv)),M=pCond,tol=p.invTol,x0=self.xo)
        #(sigma,F) = bicgstab(INTww,transpose(mat(RHSv)),M=pCond,tol=p.invTol,x0=self.xo)
        sigma = list(sigma)
        if F != 0:
            merr('Iterative matrix inverse did not converge.')
        else:
            print('--Sucessfully completed matrix-vector product: INTww*X...')
        return sigma

    def genPeriodicElements(self,x):
        p = self.parameters
        dx = p.dx
        L = p.LT
        if p.periodicBCs == True:
            N = p.Ndomains      # Number of periodic elements either side of current element
        else:
            N = 0
        n = (2*N) + 1
        xl = []
        xr = []
        for i in xrange(n):
            for xc in x:
                xl.append( xc-(dx/2.0)+(L*(i-N)) )
                xr.append( xc+(dx/2.0)+(L*(i-N)) )
        return xl, xr

    def RHSwall(vv,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,B,K,d,rhow,hw):
        # Deconstruct the input vector to fluid elements, wall position and wall acceleration
        (vff,vwp,vwv) = deconstV(vv,chebN,Nx,Nco)
        
        # First state-space equations v1_dot = v2
        RHSw = vwv
        
        # Evaluate the fourth-order spatial derivative
        d4ndx = d4dx(append(append(0,vwp),0),1,dx)
        d4ndx = d4ndx[1:(size(d4ndx)-1)]
        
        # Get the terms for the second order v2_dot = RHSwall
        RHSw2 = (-1*((K*vwp) + (B*d4ndx) + (d*vwv)))/(rhow*hw)

        RHSw = append(RHSw,RHSw2)
        return RHSw

    def d2dyCheb(self,v,chebMat):
        # Calculate the gradients in "v" using chebyChev differentiation matrix
        nn = size(v)
        nc = size(chebMat,0)
        ncols = int(nn/nc)
        R = zeros(nn)
        for i in range(ncols):
            R[(i*nc):(i*nc)+nc] = transpose(chebMat*transpose(mat(v[(i*nc):(i*nc)+nc])))
        return R

class naiveMethod():

    def __init__(self):
        self.parameters = parameters.simulation()
        self.geometry = geometry()
        # Adjust panel tolerance on FMM call (and warn if necessary)
        g = self.geometry
        p = self.parameters
        self.panTol = 0.001
        if ((1.0-g.yc[0])/p.dx) < 0.001:
            self.panTol = ((1.0-g.yc[0])/p.dx)/2.0
            print('ADJUSTING PANEL TOLERANCE FOR FMM TO ' + str(self.panTol))

    def runOdeSolver(self):
        self.generateInfluenceCoefficients()
        self.saveVariables()
        self.generateLHSandRHS()
        self.makeGeneralized()
        
        Ny = p.chebN
        Nx = p.Nx
        if p.deterministicBCs == True:
            Ny = Ny - 2

        if p.fluidOnly == False:
            N = Ny*Nx + (p.Nco-1)*2
        else:
            N = Ny*Nx        
        
        def tfun(t,y):
            yd = [0.0]*len(y)
            for i in range(len(n.RHSGeneral)):
                r = self.RHSGeneral[i]
                yd[i] = sum(r*y)
            return yd
        yinit = [0.0]*N
##        # Initialize with a spot disturbance
##        nn = int(round(Nx/4.0)*Ny + round(Ny/2.0))
##        yinit[nn] = 1e-2;yinit[nn+Ny*20] = -1e-2
##        yinit[nn-1] = 1e-2;yinit[nn-1+Ny*20] = -1e-2
        # Initialize with a sine-wave disturbance at the centreline
        NyHalf = int(round(Ny/2.0))
        for i in xrange(Nx):
            yinit[i*Ny + NyHalf] = 1e-2*np.sin(pi*g.xc[i])
            yinit[i*Ny + NyHalf-1] = 1e-2*np.sin(pi*g.xc[i])
        # Set the timeslots
        tslot = [0,p.Nsteps]
        # Call the ode45 solver
        o = octFuncs.ode()
        path = 'results/' + str(int(self.parameters.R))
        if self.parameters.fluidOnly == True:
            path = path + '/fluidOnly/'
        else:
            path = path + '/FSI/'
        o.outPath = path
        o.ode45(tfun,tslot,yinit,MaxStep=p.maxStep,InitialStep=p.maxStep)

    def runSolver(self):
        self.generateInfluenceCoefficients()
        self.saveVariables()
        self.generateLHSandRHS()
        self.makeGeneralized()
        self.callOctave()

    def makeGeneralized(self):
        print 'Making LHS and RHS Matrices into single generalised matrix: A'
        self.RHSGeneral = linalg.solve(self.LHS,self.RHS)
        print 'Done'

    def callOctave(self):
        # Save generalised matrix to A.mat
        outDict = {}
        outDict['A'] = self.RHSGeneral
        savemat('A.mat',outDict)
        # Run octave script from command line
        os.system(self.parameters.octave + ' --eval callOctaveSolvers.m')

    def generateLHSandRHS(self):

        p = self.parameters
        g = self.geometry
        f = octFuncs.finiteDifference()
    
        Nn = p.Nco + 1

        if p.fluidOnly == True:
            print 'Generating LHS and RHS matrices FOR FLUID ONLY'
        else:
            print 'Generating LHS and RHS matrices FOR FLUID+COMPLIANT WALL'

        def columnwiseMultiply(M,v):
            n = len(v)
            for i in xrange(size(M,1)/n):
                cols = range(i*n,(i+1)*n)
                M[:,cols] = M[:,cols]*v
            return M

        # A matrix for getting averages of two points
        # (for averaging node velocities to get velocities at the centre of the panels)
        avmat = zeros((p.Nco,Nn))
        for s in xrange(p.Nco):
            avmat[s,s] = 0.5
            avmat[s,s+1] = 0.5

        if p.deterministicBCs == True:
            print 'Using Deterministic No Slip Boundary Conditions'
            
            dy = g.dy[1:p.chebN-1]
            #Nf = p.Nf - 2*p.Nx
            #Ny = p.chebN - 2
            #ycF = g.ycF[1:p.chebN-1]
            
            # Extract flow-wall influence coefficients for nearest elements
            t = [i*p.chebN for i in xrange(p.Nx)] + [(i+1)*p.chebN-1 for i in xrange(p.Nx)]
            INfwWallVorts = self.ICs['INfw'][:,t]*g.dy[0]
            ITfwWallVorts = self.ICs['ITfw'][:,t]*g.dy[0]
            
            # Solve for the wall source strengths and nearest-to-wall flow elements together
            Iww = concatenate((concatenate((self.ICs['INww'],INfwWallVorts),axis=1),concatenate((self.ICs['ITww'],ITfwWallVorts),axis=1)),axis=0)
            
            # Create Right-hand-side for fluid elements 
            R = -1.0*concatenate((self.ICs['INfw'],self.ICs['ITfw']),axis=0)
            R = columnwiseMultiply(R,g.dy)
            # Set nearest-to-wall elements to zero - these are being calculated exactly to enforce the no-slip condition
            for col in t:
                R[:,col] = R[:,col]*0.0
            
            # Add compliant panel accelerations and velocities
            WallVels = zeros((p.Nx*2,Nn))
            WallVels[p.Nx+p.Nup:p.Nx+p.Nup+p.Nco,:] = avmat
            R = concatenate((R,concatenate((WallVels,zeros((p.Nx*2,Nn))),axis=0)),axis=1)           # Velocity terms
            R = concatenate((R,concatenate((zeros((p.Nx*2,Nn)),zeros((p.Nx*2,Nn))),axis=0)),axis=1) # Acceleration terms
            Mws = linalg.solve(Iww,R)
            # Get the matrix that gives wall source element strengths when multiplied with the vector of flow
            # element strengths and compliant wall velocities 
            MA = Mws[0:2*p.Nx,:]
            MAf = Mws[2*p.Nx:4*p.Nx,:]
            
        else:
            print 'Using Vorticity Injection Method to enforce BCs'
            
            # Matrix "A" - the coefficient that gives wall element strengths when
            # pre-multiplied with the vector of flow element strengths and wall
            # velocities
            # Get -1.0*INww\INfw
            MA = -1.0*linalg.solve(self.ICs['INww'],self.ICs['INfw'])
            # Multiply by fluid-element heights
            MA = columnwiseMultiply(MA,g.dy)
            # Add compliant panel accelerations and velocities
            Tmp = zeros((p.Nx*2,p.Nx*2))
            for i in xrange(p.Nx+p.Nup,p.Nx+p.Nup+p.Nco):
                Tmp[i,i] = 1.0
            #Tmp = linalg.solve(self.ICs['INww'],Tmp)    ## Check this!  Not right!!
            Tmp = Tmp[:,p.Nx+p.Nup:p.Nx+p.Nup+p.Nco]
            Tmp = dot(Tmp,avmat)
            MA = concatenate((MA,zeros((p.Nx*2,Nn)),Tmp),axis=1)
            
        # Get the number of vertical fluid elements
        Ny = p.chebN
        
        # Matrix "B" - the coefficient that gives vertical velocity at the flow
        # elements "vp" when pre-multiplied with the vector of flow element
        # strengths and wall velocities
        MB = columnwiseMultiply(self.ICs['INff'],g.dy)
        MB = concatenate((MB,zeros((p.Nf,Nn)),zeros((p.Nf,Nn))),axis=1)
        MB = MB + dot(self.ICs['INwf'],MA)

        # Matrix "C" - Right hand side of fluid transport equation
        if p.periodicBCs == False:
            MC = f.d1dx(eye(p.Nf),len(g.dy),p.dx,order=2)
            Tmp = f.d2dx(eye(p.Nf),len(g.dy),p.dx)              # Second derivative in the x-direction
        else:
            # Set up periodic finite difference matrices
            n = p.Nf + 4*p.chebN                                # Add two columns of virtual fluid elements either side of domain
            MC = f.d1dx(eye(n),len(g.dy),p.dx,order=2)
            Tmp = f.d2dx(eye(n),len(g.dy),p.dx)                 # Second derivative in the x-direction
            nrL = range(0,2*p.chebN)                            # Rows to remove (LEFT SIDE)
            nrR = range(p.Nf+2*p.chebN,p.Nf+4*p.chebN)          # Rows to remove (RIGHT SIDE)
            MC = delete(MC,nrL+nrR,axis=0)
            Tmp = delete(Tmp,nrL+nrR,axis=0)
            MC[:,p.Nf:2*p.chebN+p.Nf] = MC[:,p.Nf:2*p.chebN+p.Nf] + MC[:,nrL]
            MC[:,2*p.chebN:4*p.chebN] = MC[:,2*p.chebN:4*p.chebN] + MC[:,nrR]
            Tmp[:,p.Nf:2*p.chebN+p.Nf] = Tmp[:,p.Nf:2*p.chebN+p.Nf] + Tmp[:,nrL]
            Tmp[:,2*p.chebN:4*p.chebN] = Tmp[:,2*p.chebN:4*p.chebN] + Tmp[:,nrR]
            MC = delete(MC,nrL+nrR,axis=1)
            Tmp = delete(Tmp,nrL+nrR,axis=1)
            
        for i in xrange(p.Nf):
            s = mod(i,Ny)
            U = -1.0*(1 - (g.ycF[s]**2.0))
            MC[i,:] = MC[i,:]*U
        # Add zeros for the wall acceleration and velocity
        MC = concatenate((MC,zeros((p.Nf,Nn)),zeros((p.Nf,Nn))),axis=1)
        # Add the component due to y-direction mean-vorticity transport
        MC = MC + MB*2.0
        # Add diffusion terms
        # --- Second derivative in x-direction calculated above as "Tmp"
        # Second derivative in the y-direction - using Chebychev difference matrix
        for i in xrange(p.Nx):
            st = i*Ny
            fin = (i+1)*Ny
            Tmp[st:fin,st:fin] = Tmp[st:fin,st:fin] + g.d2dy
        Tmp = p.nu*Tmp
        Tmp = concatenate((Tmp,zeros((p.Nf,Nn)),zeros((p.Nf,Nn))),axis=1)
        # Add diffusion to RHS
        MC = MC + Tmp

        # Matrix "D" - Left hand side of fluid transport equation
        MD = concatenate((eye(p.Nf),zeros((p.Nf,Nn)),zeros((p.Nf,Nn))),axis=1)
        if p.allowSlip == False and p.deterministicBCs == False:
            # Matrix "P" - the coefficient that gives tangential velocities at the wall
            # element surfaces.
            MP = columnwiseMultiply(self.ICs['ITfw'],g.dy)
            MP = concatenate((MP,zeros((p.Nx*2,Nn)),zeros((p.Nx*2,Nn))),axis=1)
            MP = MP + dot(self.ICs['ITww'],MA)
            
            for i in xrange(p.Nx):
                print 'Injecting vorticity at Nx = ' + str(i)
                t = i*p.chebN
                b = t + p.chebN - 1
                MD[t,:] = MD[t,:] - (MP[i,:]/g.dy[0])
                MD[b,:] = MD[b,:] + (MP[p.Nx+i,:]/g.dy[p.chebN-1])
        else:
            MP = columnwiseMultiply(self.ICs['ITfw'],[0.0])
            MP = concatenate((MP,zeros((p.Nx*2,Nn)),zeros((p.Nx*2,Nn))),axis=1)
        
        # Matrix "E" - Left hand side of equation for "v1_dot = v2"
        ME = concatenate((zeros((Nn,p.Nf)),eye(Nn),zeros((Nn,Nn))),axis=1)
        
        # Matrix "F" - Right hand side of equation for "v1_dot = v2"
        MF = concatenate((zeros((Nn,p.Nf)),zeros((Nn,Nn)),eye(Nn)),axis=1)

        # Matrix "G" - Left hand side of wall equation
        # First the wall inertia terms
        MG = concatenate((zeros((Nn,p.Nf)),zeros((Nn,Nn)),(p.rhow*p.hw)*eye(Nn)),axis=1)
        # Add the pressure forcing terms to LHS of wall equation
        Ntot = p.Nf + Nn*2
        if p.deterministicBCs == True:
            PF = g.dy[0]*p.dx*MAf[p.Nx+p.Nup:p.Nx+p.Nup+p.Nco,:]
            PF = concatenate((zeros((1,Ntot)),PF),axis=0)
            PF = -1.0*cumsum(PF,axis=0)
            Pav = sum(PF,axis=0)/size(PF,0)                 # Get the average pressure across compliant wall
            PF = PF - Pav                                   # Subtract this from the nodal pressure values
            MG = MG + PF
        else:
            PF = p.dx*cumsum(MP[p.Nx+p.Nup:p.Nx+p.Nup+p.Nco-1,:],axis=0)
            MG = MG + concatenate((zeros((1,Ntot)),PF,zeros((1,Ntot))),axis=0)

        # Matrix "H" - Right hand side of wall equation
        MH = f.d4dx(eye(Nn),1,p.dx)
        MH = concatenate((zeros((Nn,p.Nf)), -1.0*(p.B*MH + p.K*eye(Nn)), -1.0*p.d*eye(Nn)), axis=1)

        # Construct the left and right hand matrices
        self.LHS = concatenate((MD,ME,MG),axis=0)
        self.RHS = concatenate((MC,MF,MH),axis=0)

        # Remove the end-nodes (which don't move) from the set of equations
        nn = [p.Nf, p.Nf+Nn-1, p.Nf+Nn, p.Nf+2*Nn-1]    # Rows/columns to delete
        for i in xrange(2):
            self.LHS = delete(self.LHS,nn,axis=i)
            self.RHS = delete(self.RHS,nn,axis=i)
        MAf = delete(MAf,nn,axis=1)

        if p.deterministicBCs == True:
            # Substitute fluid elements that are calculated deterministically
            # Firstly: extract deterministically-calculated fluid elements
            A = self.RHS[:,t]
            # Remove rows and columns that are calculated deterministically
            for i in xrange(2):
                self.LHS = delete(self.LHS,t,axis=i)
                self.RHS = delete(self.RHS,t,axis=i)
            MAf = delete(MAf,t,axis=1)
            A = delete(A,t,axis=0)
            self.RHS = self.RHS + dot(A,MAf)

        if p.fluidOnly == True:
            if p.deterministicBCs == True:
                Nf = (p.chebN - 2)*p.Nx
            else:
                Nf = p.Nf
            # Get the LHS and RHS for the fluid only
            self.LHS = self.LHS[0:Nf,0:Nf]
            self.RHS = self.RHS[0:Nf,0:Nf]

        print 'Done'

        return None

    def genPeriodicElements(self,x,y):
        p = self.parameters
        dx = p.dx
        L = p.LT
        if p.periodicBCs == True:
            N = p.Ndomains      # Number of periodic elements either side of current element
        else:
            N = 0
        n = (2*N) + 1
        xl = [x-(dx/2.0)+(L*(i-N)) for i in xrange(n)]
        xr = [x+(dx/2.0)+(L*(i-N)) for i in xrange(n)]
        return xl, xr

    def generateInfluenceCoefficients(self,forceMake=False,saveResults=True):

        p = self.parameters
        g = self.geometry        
        
        # Check if influence coefficients exist and regenerate if necessary
        remakeICs = True
        if forceMake == False:
            if os.path.exists('ICs.mat') or os.path.exists('VARS.mat'):
                invec = loadmat('VARS.mat',struct_as_record=True)
                if (int(invec['Nx'])==p.Nx) and (int(invec['Ndomains'])==p.Ndomains) and (int(invec['chebN'])==p.chebN) and (int(invec['LT'])==int(p.LT)) and (str(p.periodicBCs) in invec['periodicBCs']):
                    print 'Geometry seems to be the same as that already saved...'
                    print 'Loading existing Influence Coefficient matrices...'
                    self.ICs = loadmat('ICs.mat',struct_as_record=True)
                    print 'Done'
                    remakeICs = False
            else:
                print 'ICs.mat or VARS.mat does not exist.'
                print 'Regenerating Influence Coefficient matrices.'
        
        if remakeICs == True:
            dx = p.dx
            self.ICs = {}
            
            # Get the number of elements
            Nf = len(g.xcF)
            Nw = len(g.xcW)
            
            # Generate the list of all the elements (fluid - top wall - bottom wall)
            X = list(g.xcF) + list(g.xcW)
            Y = list(g.ycF) + list(g.ycW)
            
            # Solve the velocity at all elements due to fluid elements
            self.ICs['INff'] = zeros((Nf,Nf))
            self.ICs['ITff'] = zeros((Nf,Nf))
            self.ICs['INfw'] = zeros((Nw,Nf))
            self.ICs['ITfw'] = zeros((Nw,Nf))
            for i in xrange(Nf):
                print 'Calculating fluid-all influence for fluid element ' + str(i) + ' of ' + str(Nf)
                x = g.xcF[i]
                y = g.ycF[i]
                (xl,xr) = self.genPeriodicElements(x,y)
                (U,V) = velocitylib.vortex_element_vel(xl,[y]*len(xl),[1.0]*len(xl),xr,[y]*len(xr),[1.0]*len(xr),X,Y,threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
                self.ICs['INff'][:,i] = V[0:Nf]
                self.ICs['ITff'][:,i] = U[0:Nf]
                self.ICs['INfw'][:,i] = V[Nf:Nf+(2*Nw)]
                self.ICs['ITfw'][:,i] = U[Nf:Nf+(2*Nw)]

            # Solve the velocity at all elements due to wall elements
            self.ICs['INww'] = zeros((Nw,Nw))
            self.ICs['ITww'] = zeros((Nw,Nw))
            self.ICs['INwf'] = zeros((Nf,Nw))
            self.ICs['ITwf'] = zeros((Nf,Nw))
            for i in xrange(Nw):
                print 'Calculating wall-all influence for wall element ' + str(i) + ' of ' + str(Nw)
                x = g.xcW[i]
                y = g.ycW[i]
                (xl,xr) = self.genPeriodicElements(x,y)
                (U,V) = velocitylib.source_element_vel(xl,[y]*len(xl),[1.0]*len(xl),xr,[y]*len(xr),[1.0]*len(xr),X,Y,threads=p.Nthreads,fmm=True,assume_point_length=p.apl,panel_tolerance=self.panTol)
                self.ICs['INwf'][:,i] = V[0:Nf]
                self.ICs['ITwf'][:,i] = U[0:Nf]
                self.ICs['INww'][:,i] = V[Nf:Nf+(2*Nw)]
                self.ICs['ITww'][:,i] = U[Nf:Nf+(2*Nw)]
                # Correct the self-inflence (diagonals) for wall-wall influence coefficients
                if y > 0:
                    self.ICs['INww'][i,i] = -0.5

            if saveResults == True:
                self.saveInfluenceCoefficients()
                self.saveVariables()
        return None

    def saveInfluenceCoefficients(self):
        # Save the variables to a .mat file (matlab/octave format)
        print 'Writing influence coefficient matrices to ICs.mat'
        savemat('ICs.mat',self.ICs)
        print 'Finished...'
        
    def saveVariables(self):
        import inspect
        outvec = {}
        print 'Writing variables to VARS.mat'
        vars = inspect.getmembers(self.parameters)
        for v in vars:
            if type(v[1]) == float or type(v[1]) == int:
                outvec[v[0]] = array([[v[1]]])
            elif type(v[1]) == str:
                outvec[v[0]] = v[1]
            elif type(v[1]) == bool:
                outvec[v[0]] = str(v[1])
        savemat('VARS.mat',outvec)
        print 'Finished'

        return None





class postProc():

    def __init__(self):
        self.ps = parameters.simulation()
        self.p = parameters.postProcessing()
        self.g = geometry()

    def incLineType(self,k):
        if k == 0:
            lt = 'k-'
        elif k == 1:
            lt = 'k--'
        elif k == 2:
            lt = 'k-.'
        elif k == 3:
            lt = 'k:'
        elif k == 4:
            lt = 'k-o'
        elif k == 5:
            lt = 'k-v'
        else:
            lt = 'k-'
        return lt

    def getResultsDir(self,R = []):
        if R == []:
            R = int(self.ps.R)
        Rpath = 'results/' + str(R)
        if self.ps.fluidOnly == True:
            Rpath = Rpath + '/fluidOnly/'
        else:
            Rpath = Rpath + '/FSI/'
        return Rpath

    def getFilesByPrefix(self,prefix,Rpath):
        allFiles = os.listdir(Rpath)
        files = []
        for f in allFiles:
            if prefix in f:
                files.append(f)
        files.sort()
        return files

    def plotGridY(self):
        
        fig1 = plt.figure()
        plt.hold(True)
        plt.plot([0.0]*len(g.yn),g.yn,'kx-')
        plt.plot([0.0]*len(g.yc),g.yc,'go')
        plt.hold(False)
        
        plt.show()

class ppOde45(postProc):
    
    def plotMonitors(self,eleList,RList=[],ax=[],leg=[],legLoc=0,saveMats=False):
        '''
        RList is a list of integers or strings that state the directory in which to
        extract results from within the "results" directory.
        
        eleList is a list of integers or a list of lists of integers which are the elements
        from which to extract the useful information.  NOTE: if RList is undefined then
        len(eleList) = 1 else the length of eleList must equal the length of Rlist.
        '''        

        if RList == []:
            RList = [int(self.ps.R)]

        fig1 = plt.figure()
        plt.xlabel('Time (s)')
        plt.ylabel('Perterbation vorticity at monitor point')
        a = plt.axes()
        a.ticklabel_format(axis='y',scilimits=(-1,1))
        plt.hold(True)
        
        M = [[] for i in xrange(len(RList))]

        k = 0
        for R in RList:
            print 'Getting monitor points for R = ' + str(R)
            Rpath = self.getResultsDir(R)
            fnames = self.getFilesByPrefix('TStep',Rpath)
            f = 0

            if type(eleList[k]) == list:
                n = len(eleList[k])
            else:
                n = 1
                eleList[k] = [eleList[k]]

            Z = [[] for i in xrange(n)]
            t = []
            for fname in fnames:
                print 'Loading frame: ' + fname
                fpath = Rpath + fname
                inDict = loadmat(fpath,struct_as_record=True)
                v = inDict['y']
                tt = inDict['t']
                
                t.append(tt[0][0])
                for i in xrange(n):
                    Z[i].append(v[eleList[k][i]][0])
                
            M[k] = Z
            for i in xrange(n):
                plt.plot(t,Z[i],self.incLineType(k))
            k += 1

        plt.hold(False)
        plt.grid(True)
        
        if ax != []:
            plt.axis(ax)
        if leg != []:
            plt.legend(leg,loc=legLoc)
        fig1.savefig('plotMonitors.eps')

        if saveMats == True:
            print 'Saving matrix to file monitors.mat'
            o = {}
            o['M'] = M
            savemat('monitors.mat',o)
    
    def makeMovie(self,prt='Fluid',R=[],fillConts=False,Nlevels=20,dumpFigs=True,keepFigs=False):
        '''
        prt = 'Fluid' or 'Wall' to plot wall or fluid displacements
        fillConts - plot with contourf or contour
        '''

        if R == []:
            R = int(self.ps.R)

        p = self.p.movie()
        g = self.g
        pg = self.p.general()
        fig1 = plt.figure()
        fig1.set_figwidth(16);fig1.set_figheight(4)

        Ny = self.ps.chebN
        yc = g.yc
        if self.ps.deterministicBCs == True:
            yc = g.yc[1:Ny-1]
            Ny = self.ps.chebN - 2
        
        Nx = self.ps.Nx
        if self.ps.fluidOnly == True:
            Ncw = 0
        else:
            Ncw = (self.ps.Nco + 1)*2

        Nco = self.ps.Nco
        Rpath = self.getResultsDir(R)
        fnames = self.getFilesByPrefix('TStep',Rpath)
        f = 0
        Z = []
        alim = 1.0e-8
        for fname in fnames:
            fpath = Rpath + fname
            inDict = loadmat(fpath,struct_as_record=True)
            if prt == 'Fluid':
                vm = inDict['y'][0:Ny*Nx]
                vmT = reshape(vm,(Ny,Nx),order='F')
                
                m = np.max(real(vmT))
                m = m*pg.maxFactor;
                if Z == []: 
                    Z = np.linspace(-1.0*m,m,Nlevels)
                elif m > 2*np.max(Z) or m < 0.5*np.max(Z):
                    Z = np.linspace(-1.0*m,m,Nlevels)

                plt.xlabel('x-position (m)')
                plt.ylabel('y-position (m)')
                if fillConts == True:
                    plt.contourf(g.xc,yc,vmT,Z)
                    plt.axis('tight')
                    plt.colorbar()
                else:
                    l = ['dashed']*(Nlevels/2) + ['solid']*(Nlevels/2)
                    plt.contour(g.xc,yc,vmT,Z,colors='k',linestyles=l)
                    plt.axis('tight')
                    plt.grid(True)
            else:
                vm = inDict['y'][Ny*Nx:Ny*Nx+Nco-1]
                vm = array([0.0] + [vv for vv in vm] + [0.0])
                plt.plot(vm,'k-')
                absmax = max(abs(vm))
                if absmax > alim:
                    alim = absmax
                plt.axis([0,len(vm)-1,-1.0*alim,alim])
                plt.grid(True)

            if dumpFigs == True:
                fn = '_tmp%04d.png'%f
                print 'Saving frame', fn
                fig1.savefig(fn)
            fig1.clf()
            f += 1
        
        print 'Making movie animation.mpg - this make take a while'
        os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(int(p.fps)) + " -ovc lavc -lavcopts vcodec=" + p.vcodec + " -nosound -o ode45_" + prt + ".avi")
        if keepFigs == False:        
            os.system("rm -v *.png")

class ppEigs(postProc):

    def makeMovie(self,fname,prt='Fluid',modeNum=0,fillConts=False,Nlevels=20,dumpFigs=True,keepFigs=False):
        
        pg = self.p.general()
        pm = self.p.movie()
        ps = parameters.simulation()

        Np = pm.Np
        P = pm.P

        x = self.g.xc
        if self.ps.deterministicBCs == True:
            y = self.g.yc[1:len(self.g.yc)-1]
        else:
            y = self.g.yc
        x = array(x)
        y = array(y)
        L = ps.LT
        Nf = len(y)*len(x)
        Nw = (ps.Nco-1)*2

        # NOTE: figsize gives the figure size in inches @ 100dpi => 16=1200pixels
        fig1 = plt.figure(figsize=(16,4))
        #ax = fig1.add_subplot(111)

        print fname
        indict = loadmat(pg.resultsDir + '/' + fname,struct_as_record=True)
        v = [indict['Veigs'][i][modeNum] for i in xrange(len(indict['Veigs']))]
        e = indict['evals'][modeNum]
        if pg.ignoreTemporalGrowth == True:
            E = complex(0,imag(e))
        else:
            E = e
        print 'Eigenvalue = ' + str(E)

        if prt == 'Fluid':
            vm = v[0:Nf]
            vm = reshape(vm,(len(y),len(x)),order='F')
        else:
            vm = array(v[Nf:Nf+Nw])
            vm = vm[0:Nw/2]
            vm = concatenate(([0.0],vm,[0.0]))

        T = abs(2*pi/imag(e))
        tRange = linspace(0,T*P,Np*P,endpoint=pm.includeEndTime)
        Mv = [[] for i in xrange(len(tRange))]
        f = 0
        Z = []
        alim=1e-8
        for t in tRange:
            timeComponent = exp(E*t)
            vmT = vm*timeComponent
            vmT = real(vmT)
            if prt == 'Fluid':
                if pg.normalizeSpatialGrowth == True:
                    relativeSpatialGrowth = zeros((len(y),len(x)))
                    for i in xrange(len(y)):
                        U = 1 - (y[i]**2)
                        #U = 1.0
                        relativeSpatialGrowth[i,:] = exp(real(e)*(x/U))
                    vmT = multiply(vmT,relativeSpatialGrowth)

                m = np.max(real(vmT))
                m = m*pg.maxFactor;
                if Z == []: 
                    Z = np.linspace(-1.0*m,m,20)
                elif m > 2*np.max(Z) or m < 0.5*np.max(Z):
                    Z = np.linspace(-1.0*m,m,20)

                if fillConts == True:
                    plt.contourf(x,y,vmT,Z)
                    plt.axis('tight')
                    plt.colorbar()
                else:
                    l = ['dashed']*(Nlevels/2) + ['solid']*(Nlevels/2)
                    plt.contour(x,y,vmT,Z,colors='k',linestyles=l)
                    plt.axis('tight')
                    plt.grid(True)
            else:
                plt.plot(vmT,'k-')
                absmax = max(abs(vmT))
                if absmax > alim:
                    alim = absmax
                plt.axis([0,len(vmT)-1,-1.0*alim,alim])
                plt.grid(True)
                Mv[f] = vmT
            
            if dumpFigs == True:
                fn = '_tmp%04d.png'%f
                print 'Saving frame', fn
                fig1.savefig(fn)
            f += 1
            fig1.clf()
        
        if prt != 'Fluid':
            a = plt.axes()
            a.ticklabel_format(axis='y',scilimits=(-1,1))
            plt.hold(True)
            for i in xrange(len(Mv)):
                plt.plot(Mv[i],'k-')
            plt.hold(False)
            plt.axis([0,len(vmT)-1,-1.0*alim,alim])
            plt.grid(True)
            fig1.savefig('wallSnapshots.eps')

        print 'Making movie animation.mpg - this make take a while'
        os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(int(pm.fps)) + " -ovc lavc -lavcopts vcodec=" + pm.vcodec + " -nosound -o eigs_" + prt + ".avi")
        if keepFigs == False:        
            os.system("rm -v *.png")


if __name__ == "__main__":
    p = parameters.simulation()
    if p.method == 'Naive':
        print 'Using NAIVE method in Python to construct square matrices...'
        n = naiveMethod()
        if p.solver == 'eigs':
            print 'Solving EIGENVALUE problem...'
            n.runEigSolver()
        elif p.solver == 'ode45':
            print 'Solving Intial Value problem (ODE45)...'
            n.runOdeSolver()
        else:
            merr('Solver type not recognized.  Must be eigs or ode45.')
    elif p.method == 'FMM':
        print 'Using FMM method...'
        f = fmmMethod()
        if p.solver == 'eigs':
            print('Solving EIGENVALUE problem...')
            f.runEigSolver()
        elif p.solver == 'ode45':
            print 'Solving TIME-STEPPING (inital value) problem...'
            f.runOdeSolver()
##            # Profiler (for optimisation)
##            prof = cProfile.Profile()
##            prof = prof.runctx("f.runOdeSolver()", globals(), locals())
##            stats = pstats.Stats(prof)
##            stats.sort_stats("cumulative")  # Or cumulative
##            stats.print_stats(80)  # 80 = how many to print
##            # ----------------
##        if os.path.exists('v.mat'):
##            # Variable thrown from octave already exists so go straight to multiplication
##            inDict = loadmat('v.mat')
##            v = inDict['v']
##            n = fmmMethod()
##            x = n.RHSfluid(v)
##            outDict = {}
##            outDict['x'] = x
##            savemat('x.mat',outDict)
    else:
        merr('Method type not recognized.  Must be Naive or FMM')


    

