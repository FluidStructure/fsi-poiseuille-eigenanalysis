from scipy import reshape, real, imag, array, exp, pi, linspace, divide, multiply
from scipy import append, size, concatenate, dot, eye, mod, linalg
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib.pyplot import imshow, axis, contour, contourf
from numpy import mat, reshape, linspace, reshape, append, random, ones, delete
from numpy import meshgrid, bmat, array, size, zeros, transpose, cumsum
import numpy as np
import os, sys, subprocess

import pdb

import sys
sys.path.append('../python_tools.git/')
import velocitylib

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
            self.method = p[14]
            self.octave = p[15]
            self.solver = p[16]
            
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

        c = chebychev()

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
    
    def runFMM(vv,sigma):
        # Read in parameters for the system
        (nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,B,K,d,rhow,hw,dy) = getParams()

        def multRHS(x):
            if iscomplexobj(x) == False:
                #print 'x is only REAL!'
                #print 'CALCULATING RHS*x!'
                RHSf = RHSfluid(x,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,dy)
            else:
                #print 'x is COMPLEX!'
                #print 'CALCULATING RHS*real(x)!'
                RHSfr = RHSfluid(real(x),nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,dy)
                #print 'CALCULATING RHS*imag(x)!'
                RHSfc = RHSfluid(imag(x),nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,dy)
                RHSf = RHSfr + RHSfc
            return RHSf
        RHS = LinearOperator( (Nx*(chebN-2),Nx*(chebN-2)), matvec=multRHS, dtype='float64' )
        # Calculate the RHS*x of the fluid equations
        # Same as: (RHSf,uSlip) = RHSfluid(vv,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,dy)

        if sigma == 0.0:
            xoutRe = RHS.matvec( vv )
            xoutIm = 0*xoutRe
        else:
            # OR, RHS\x (when sigma is specified ~= 0.0 or sigma = 'SM'):
            #(xout,F) = gmres(RHS, vv, tol=1.0000000000000001e-02, restrt=50, maxiter=100)
            #(xout,F) = cgs(RHS, vv, tol=1.0000000000000001e-02, maxiter=100)    
            print 'sigma = ', sigma
            F = 1;
            Ny = chebN-2;
            #pdb.set_trace()
            v = append(real(vv),imag(vv))
            
            def preCondRHS(xx):
                xo = xx
                return xo
            pCondRHS = LinearOperator( (2*Nx*Ny,2*Nx*Ny), matvec=preCondRHS, dtype='float64' )
            def multRHScomplex(xi):
                #pdb.set_trace()
                xiR = xi[0:Nx*Ny];
                xiC = xi[Nx*Ny:2*Nx*Ny]
                xoR = (multRHS(xiR) + real(-1*sigma)*xiR) - imag(-1*sigma)*xiC
                xoC = imag(-1*sigma)*xiR + (multRHS(xiC) + real(-1*sigma)*xiC)
                xo = append(xoR,xoC)
                return xo
            RHScomp = LinearOperator( (2*Nx*Ny,2*Nx*Ny), matvec=multRHScomplex, dtype='float64' )
            #(xout,F) = minres(RHScomp, v, M=None, tol=1.0e-6, maxiter=1000,shift=complex(real(sigma),imag(sigma))) 
            (xo,F) = minres(RHScomp, v, M=pCondRHS, tol=1.0e-3, maxiter=1000)
            xoutRe = xo[0:Nx*Ny]
            xoutIm = xo[Nx*Ny:2*Nx*Ny] 
            if F == 0:
                print 'SUCESSFULLY COMPLETED RHScomplex\X!'
            else:
                pdb.set_trace()

        # Calculate the RHS of the wall equations
        # RHSw = RHSwall(vv,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,B,K,d,rhow,hw)
        
        # Combine the RHSwall and RHSfluid results
        # RHS = append(RHSf,RHSw)
        
        # Get the inv(LHS)*(RHS*x)
        # x = invLHS(RHS,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,rhow,hw)
        
        return xoutRe,xoutIm,pcyy[0:chebN-2]

    def invLHS(RHS,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,rhow,hw):

        def dummyLHS(vl):
            print 'Inside dummyLHS...'
            X = vl
            for i in range(Nx*chebN+(Nco-1),Nx*chebN+2*(Nco-1)):
                X[i] = X[i]*5.0
            return X

        def multLHS(vl):
            #pdb.set_trace()
            # Evaluating the matrix-vector product with the LHS matrix of wall and fluid equations with input (x)
            print 'Evaluating matrix-vector product: LHS*vf...'

            # Deconstruct the input vector to fluid elements, wall position and wall acceleration
            (vff,vwp,vwv) = deconstV(vl,chebN,Nx,Nco)

            # Get the velocity of each wall panel
            wpv = panelVels(vwv,Nup,Nco,Ndn)

            # Solve the velocity flux through wall elements due to fluid elements
            (fwu,fwv) = velocitylib.vortex_element_vel(pcxx-(dx/2),pcyy,vff,pcxx+(dx/2),pcyy,vff,xw,yw)
            
            # Solve the wall element strengths
            def multINww(x):
                # Evaulate the matrix-vector product of wall-normal influence coefficients using Jarrads FMM
                print 'Evaluating matrix-vector product: INww*Gamma...'
                (u,v) = velocitylib.source_element_vel(xw-(dx/2),yw,x,xw+(dx/2),yw,x,xw,yw)
                return v
            INww = LinearOperator( (Nx*2,Nx*2), matvec=multINww, dtype='float64' )
            (sigma1,Fi) = minres(INww,wpv)
            if Fi == 0:
                print 'Sucessfully completed inv(INww)*wpv'
            (sigma2,Fi) = minres(INww,fwv)
            if Fi == 0:
                print 'Sucessfully completed inv(INww)*fwv'
            sigma = sigma1 + sigma2
            
            # Evaluate the perturbation at all wall elements due to wall elements
            (wwu,wwv) = velocitylib.source_element_vel(xw-(dx/2),yw,sigma,xw+(dx/2),yw,sigma,xw,yw)
            
            # Get total slip velocity at the wall
            us = array(wwu) + array(fwu)
            
            # Initialize the return vector
            nn = chebN*Nx
            #R = array([0.0]*size(vv))      # For debugging.  Really R = vv (initialise with ones on diagonal)
            R = vl
            
            # Apply vorticity injection terms
            if False:
                for i in range(0,Nx):
                    R[i*chebN + 0] = R[i*chebN + 0] - us[i]
                    R[i*chebN + chebN-1] = R[i*chebN + chebN-1] + us[i+Nx]
            
            # Add pressure integrals from leading edge (to compliant panels)
            if True:
                pn = array([0.0]*2)
                pav = (pn[0] + pn[1])/2
                co = (Nco-1)-1
                for i in range(Nx,Nx+Nup+Nco):
                    pn[0] = pn[1]
                    pn[1] = pn[0] + us[i]*dx                  # Integrate the pressure across nodes
                    if i >= Nx+Nup:
                        pav = (pn[0] + pn[1])/2               # Average pressure across panel
                        R[Nx*chebN+co] = R[Nx*chebN+co] + (pav*dx)/(rhow*hw)
                        co = co + 1
            
            # pdb.set_trace()
            return R
        
        #  Evaluate 'x' = inv(LHS)*(RHS*v)
        #pdb.set_trace()
        nn = size(RHS)
        LHS = LinearOperator( (nn,nn), matvec=multLHS, dtype='float64' )
        #LHS = LinearOperator ( (nn,nn), matvec=dummyLHS, dtype='float64' )
        (X,F) = minres(LHS,RHS,x0=RHS)
        #pdb.set_trace()
        if F == 0:
            print 'Sucessfully completed inv(LHS)*RHS'
        
        return X

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

    def RHSfluid(vv,nu,dx,Nx,Nup,Nco,Ndn,chebN,d2dy,xw,yw,pcxx,pcyy,dy):
        
        Nthreads = 2
        fmmOn = True
        
        # Get the position of wall-vorticies (the vortex elements closest to the wall)
        xwv = xw
        ywv = (-1.0*yw*(dy[0]/2.0)) + yw
        
        # Deconstruct the input vector to fluid elements, wall position and wall acceleration
        Ny = chebN - 2;dyf = dy[1:(size(dy)-1)]
        vff = array([0.0]*size(vv))
        for i in range(0,Nx):
            vff[i*Ny:(i+1)*Ny] = vv[i*Ny:(i+1)*Ny]*array(dyf)

        # Solve the velocity flux through all elements due to fluid elements (except near-wall fluid elements)
        (fau,fav) = velocitylib.vortex_element_vel(pcxx-(dx/2),pcyy,vff,pcxx+(dx/2),pcyy,vff,append(pcxx,xw),append(pcyy,yw),threads=Nthreads,fmm=fmmOn)
        ffv = fav[0:(Ny*Nx)]
        fwv = fav[(Ny*Nx):(Ny*Nx + Nx*2)]
        fwu = fau[(Ny*Nx):(Ny*Nx + Nx*2)]

        # Solve the wall element (source and nearest vorticity) strengths (using a pre-conditioner)
        def preCond(xx):
            xo = xx
            xo[0*Nx:1*Nx] = (xo[0*Nx:1*Nx])*-2.0                # Upper source elements
            xo[1*Nx:2*Nx] = (xo[1*Nx:2*Nx])*2.0                 # Lower source elements
            xo[2*Nx:3*Nx] = (xo[2*Nx:3*Nx])*-2.0*(1/dy[0])      # Upper vortex elements
            xo[3*Nx:4*Nx] = (xo[3*Nx:4*Nx])*2.0*(1/dy[0])       # Lower vortex elements
            return xo
        def multINTww(xx):
            # Evaulate the matrix-vector product of source/sink influence coefficients using Jarrads FMM
            #print 'Evaluating matrix-vector product: INTww*X...'
            ss = xx[0:(Nx*2)]
            (us,vs) = velocitylib.source_element_vel(xw-(dx/2),yw,ss,xw+(dx/2),yw,ss,xw,yw,threads=Nthreads,fmm=fmmOn)
            # Adjust the upper panels to have a self-influence of -0.5
            # pdb.set_trace()
            #vs[0:Nx] = subtract(vs[0:Nx],ss[0:Nx])
            for i in range(0,Nx):
                vs[i] = vs[i] - ss[i]
            
            # Evaluate the matrix-vector product of wall-vortex influence coefficients
            sv = xx[(Nx*2):(Nx*4)]*dy[0]
            (uv,vv) = velocitylib.vortex_element_vel(xwv-(dx/2),ywv,sv,xwv+(dx/2),ywv,sv,xw,yw,threads=Nthreads,fmm=fmmOn)
            # Construct the output vector
            ov = append(add(vs,vv),add(us,uv))
            return ov
        INTww = LinearOperator( (Nx*2*2,Nx*2*2), matvec=multINTww, dtype='float64' )
        pCond = LinearOperator( (Nx*2*2,Nx*2*2), matvec=preCond, dtype='float64' )
        RHSw = append(multiply(fwv,-1.0),multiply(fwu,-1.0))
        (sigma,F) = minres(INTww,transpose(mat(RHSw)),M=pCond)
        if F != 0:
            pdb.set_trace()
            #print 'Sucessfully completed inv(INTww)*[fwv;fwu]'
        # pdb.set_trace()       # To test result try: INww.matvec(sigma) + fwv

        # Evaluate the perturbation velocity in all elements from the wall sources and vortices
        ss = sigma[0*Nx:2*Nx]
        sv = sigma[2*Nx:4*Nx]*dy[0]
        (sau,sav) = velocitylib.source_element_vel(xw-(dx/2),yw,ss,xw+(dx/2),yw,ss,append(pcxx,xw),append(pcyy,yw),threads=Nthreads,fmm=fmmOn)
        (vau,vav) = velocitylib.vortex_element_vel(xwv-(dx/2),ywv,sv,xwv+(dx/2),ywv,sv,append(pcxx,xw),append(pcyy,yw),threads=Nthreads,fmm=fmmOn)
        wfv = array(sav[0:(Ny*Nx)]) + array(vav[0:(Ny*Nx)])

        # Get the mean velocity profile and gradients
        # Plane Pousille flow
        U = (1 - pcyy**2)
        d2Udy2 = -2.0        # This is actually the gradient (in y-direction) of VORTICITY - don't be fooled by variable name

        # Get the vorticity-gradient multiplied with the y-direction perturbation velocity
        vpdody = (wfv + array(ffv))*d2Udy2
        
        # Get the mean velocity (1m/s in this case) multiplied with the x-direction gradient of the perturbation vorticity
        Udodx = U*array(d1dx(vv,Ny,dx,2))

        # Constuct a new vector of fluid element strengths that includes the near-wall fluid elements (determined)
        sv = sigma[2*Nx:4*Nx]
        vvw = [0.0]*chebN*Nx
        for i in range(0,Nx):
            vvw[i*chebN:(i+1)*chebN] = append(append(array(sv[i]),vv[i*Ny:(i+1)*Ny]),array(sv[i+Nx]))

        # Get the second derivatives in x and y direction multiplied by kinematic viscosity
        d2odyw = d2dyCheb(vvw,d2dy)
        d2odx = d2dx(vv,Ny,dx)

        # Extract relevant parts of d2odyw
        d2ody = [0.0]*Ny*Nx
        for i in range(0,Nx):
            d2ody[i*Ny:(i+1)*Ny] = d2odyw[(i*chebN + 1):((i+1)*chebN - 1)]
        
        # Evaluate the RHS of the fluid equation
        RHSf = nu*(array(d2ody) + array(d2odx)) - (vpdody + array(Udodx))
        #RHSf = nu*(array(d2ody) + array(d2odx))
        
        return RHSf

    def d2dyCheb(v,chebMat):
        # Calculate the gradients in "v" using chebyChev differentiation matrix
        nn = size(v)
        nc = size(chebMat,0)
        ncols = int(nn/nc)
        R = zeros(nn)
        for i in range(0,ncols):
            R[(i*nc):(i*nc)+nc] = transpose(chebMat*transpose(mat(v[(i*nc):(i*nc)+nc])))
        return R

    def deconstV(vv,chebN,Nx,Nco):
        vff = vv[0:(chebN*Nx)]
        vwp = vv[(chebN*Nx):(chebN*Nx+Nco-1)]
        vwv = vv[(chebN*Nx+Nco-1):(chebN*Nx+2*Nco-2)]
        return vff,vwp,vwv

    def panelVels(vwv,Nup,Nco,Ndn):
        Nx = Nup+Nco+Ndn
        wpv = [0.0]*Nx                                              # Upper wall panels
        wpv = append(wpv,[0.0]*Nup)                                 # Panels upstream of the compliant section
        ve = append(append(0.0,vwv),0.0)
        wpv = append(wpv,(ve[0:(size(ve)-1)] + ve[1:size(ve)])/2)   # Compliant wall section
        wpv = append(wpv,[0.0]*Ndn)                                 # Panels downstream of compliant section
        return wpv


class naiveMethod():

    def __init__(self):
        self.parameters = parameters.simulation()
        self.geometry = geometry()

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
        f = finiteDifference()
    
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
            R = columnwiseMultiply(R,dy)
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
            
            ## Remove nearest-to-wall flow elements from other influence matrices
            #self.ICs['INff'] = delete(self.ICs['INff'],t,axis=0)
            #self.ICs['INff'] = delete(self.ICs['INff'],t,axis=1)
            #self.ICs['INwf'] = delete(self.ICs['INwf'],t,axis=0)
           
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
        MC = f.d1dx(eye(p.Nf),len(g.dy),p.dx,order=2)
        for i in xrange(p.Nf):
            s = mod(i,Ny)
            U = -1.0*(1 - (g.ycF[s]**2.0))
            MC[i,:] = MC[i,:]*U
        # Add zeros for the wall acceleration and velocity
        MC = concatenate((MC,zeros((p.Nf,Nn)),zeros((p.Nf,Nn))),axis=1)
        # Add the component due to y-direction mean-vorticity transport
        MC = MC - MB*2.0
        # Add diffusion terms
        Tmp = f.d2dx(eye(p.Nf),len(g.dy),p.dx)              # Second derivative in the x-direction
        # Second derivative in the y-direction
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
        # Add the pressure forcing terms
        PF = p.dx*cumsum(MP[p.Nx+p.Nup:p.Nx+p.Nup+p.Nco-1,:],axis=0)
        Ntot = p.Nf + Nn*2
        MG = MG + concatenate((zeros((1,Ntot)),PF,zeros((1,Ntot))),axis=0)

        # Matrix "H" - Right hand side of wall equation
        MH = f.d4dx(eye(Nn),1,p.dx)
        MH = concatenate((zeros((Nn,p.Nf)), -1.0*(p.B*MH + p.K*eye(Nn)), -1.0*p.d*eye(Nn)), axis=1)

        # Construct the left and right hand matrices
        self.LHS = concatenate((MD,ME,MG),axis=0)
        self.RHS = concatenate((MC,MF,MH),axis=0)

        # Remove the end-nodes (which don't move) from the set of equations
        nn = [p.Nf, p.Nf+Nn-1, p.Nf+Nn, p.Nf+2*Nn]    # Rows/columns to delete
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

    def generateInfluenceCoefficients(self,forceMake=False,saveResults=True):

        p = self.parameters
        g = self.geometry        
        
        # Check if influence coefficients exist and regenerate if necessary
        remakeICs = True
        if forceMake == False:
            if os.path.exists('ICs.mat') or os.path.exists('VARS.mat'):
                invec = loadmat('VARS.mat',struct_as_record=True)
                if (int(invec['Nx'])==p.Nx) and (int(invec['chebN'])==p.chebN) and (int(invec['LT'])==int(p.LT)):
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
                (U,V) = velocitylib.vortex_element_vel([x-(dx/2.0)],[y],[1.0],[x+(dx/2.0)],[y],[1.0],X,Y)
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
                (U,V) = velocitylib.source_element_vel([x-(dx/2.0)],[y],[1.0],[x+(dx/2.0)],[y],[1.0],X,Y)
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
        i = 0;R[i] = 0.0 - 4.0*0.0 + 0.0*n[i] - 4.0*n[i+N] + n[i+2*N]
        i = 1;R[i] = 0.0 - 4.0*n[i-N] + 5.0*n[i] - 4.0*n[i+N] + n[i+2*N]
        i = nn-2;R[i] = n[i-2*N] - 4.0*n[i-N] + 5.0*n[i] - 4.0*n[i+N] + 0.0
        i = nn-1;R[i] = n[i-2*N] - 4.0*n[i-N] + 0.0*n[i] - 0.0 + 0.0
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


class postProc():

    def __init__(self):
        self.p = parameters.postProcessing()
        self.g = geometry()

    def getFilesByPrefix(self,prefix):
        allFiles = os.listdir('results')
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
    
    def makeMovie(self,dumpFigs=True):
        p = self.p.movie()
        g = self.g
        fig1 = plt.figure()
        fig1.set_figwidth(16);fig1.set_figheight(4)
        vars = loadmat('VARS.mat',struct_as_record=True)

        Ny = int(vars['chebN'])
        yc = g.yc
        if 'True' in vars['deterministicBCs']:
            yc = g.yc[1:Ny-1]
            Ny = int(vars['chebN']) - 2
            
        Nx = int(vars['Nx'])
        if 'True' in str(vars['fluidOnly']):
            Ncw = 0
        else:
            Ncw = (int(vars['Nco']) + 1)*2
        
        fnames = self.getFilesByPrefix('TStep')
        f = 0
        for fname in fnames:
            fpath = 'results/'+fname
            inDict = loadmat(fpath,struct_as_record=True)
            vmT = reshape(inDict['y'],(Ny,Nx),order='F')
            
            plt.contourf(g.xc,yc,vmT)
            plt.axis('tight')
            plt.colorbar()
              
            if dumpFigs == True:
                fn = '_tmp%03d.png'%f
                print 'Saving frame', fn
                fig1.savefig(fn)
            fig1.clf()
            f += 1
        
        print 'Making movie animation.mpg - this make take a while'
        os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(int(p.fps)) + " -ovc lavc -lavcopts vcodec=" + p.vcodec + " -nosound -o ode45.avi")
        #os.system("rm -v *.png")

class ppEigs(postProc):

    def plotFrameZero(self):
        #for fname in os.listdir(resultsDir):
        for fname in ['evals_R5000.mat']:
            self.plotTimes(fname,1,1)
        show()
    
    def makeMovie(self):
        p = self.p.movie()
        #for fname in os.listdir(resultsDir):
        for fname in ['evals_R5000.mat']:
            self.plotTimes(fname,p.Np,p.P,dumpFigs=True)
        print 'Making movie animation.mpg - this make take a while'
        os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(int(self.p.movie.fps)) + " -ovc lavc -lavcopts vcodec=" + self.p.movie.vcodec + " -nosound -o eigs.avi")
        #os.system("rm -v *.png")

    def plotTimes(self,fname,Np,P,dumpFigs=False):

        vdict = loadmat('VARS.mat',struct_as_record=True)
        x = [vdict['pcx'][0][i] for i in xrange(len(vdict['pcx'][0]))]
        y = [vdict['pcy'][i][0] for i in xrange(len(vdict['pcy']))]
        x = array(x)
        y = array(y)
        L = float(vdict['LT'][0][0])

        # NOTE: figsize gives the figure size in inches @ 100dpi => 16=1200pixels
        fig = plt.figure(figsize=(16,4))
        ax = fig.add_subplot(111)

        print fname
        indict = loadmat(resultsDir + '/' + fname,struct_as_record=True)
        v = [indict['Veigs'][i][modeNum] for i in xrange(len(indict['Veigs']))]
        e = indict['evals'][modeNum][0]
        if ignoreTemporalGrowth == True:
            E = complex(0,imag(e))
        else:
            E = e
        print 'Eigenvalue = ' + str(e)

        vm = reshape(v,(len(y),len(x)),order='F')

        T = abs(2*pi/imag(e))
        tRange = linspace(0,T*P,Np*P,endpoint=includeEndTime)
        f = 0
        Z = []
        for t in tRange:
            ax.cla()
            
            timeComponent = exp(E*t)
            vmT = vm*timeComponent
            vmT = real(vmT)
            if normalizeSpatialGrowth == True:
                relativeSpatialGrowth = exp(real(e)*x)
                vmT = multiply(vmT,relativeSpatialGrowth)

            if Z == []:
                m = np.max(real(vmT))
                m = m*maxFactor;
                Z = np.linspace(-1.0*m,m,20)

            plt.contourf(x,y,vmT,Z)
            plt.axis('tight')
            
            if dumpFigs == True:
                fn = '_tmp%03d.png'%f
                print 'Saving frame', fn
                fig.savefig(fn)
            f += 1


if __name__ == "__main__":
    p = parameters.simulation()
    if p.method == 'Naive':
        print 'Using NAIVE method in Python to construct square matrices...'
        n = naiveMethod()
        if p.solver == 'eigs':
            print 'Solving EIGENVALUE problem...'
            n.runSolver()
        elif p.solver == 'ode45':
            print 'Solving Intial Value problem (ODE45)...'
            n.runSolver()
        else:
            merr('Solver type not recognized.  Must be eigs or ode45.')
    elif p.method == 'FMM':
        if os.path.exists('v.mat'):
            # Variable thrown from octave already exists so go straight to multiplication
            inDict = loadmat('v.mat')
            v = inDict['v']
            n = fmmMethod()
            x = n.RHSfluid(v)
            outDict = {}
            outDict['x'] = x
            savemat('x.mat',outDict)
        else:
            # Start a new simulation from scratch
            print 'Solving TIME-STEPPING (inital value) problem...'
            print 'Using FMM method...'
    else:
        merr('Method type not recognized.  Must be Naive or FMM')


    

