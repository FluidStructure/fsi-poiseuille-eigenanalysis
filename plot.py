from scipy import reshape, real, imag, array, exp, pi, linspace, divide, multiply
from scipy.io import loadmat
from matplotlib.pyplot import figure, show
from matplotlib.image import NonUniformImage
from matplotlib.pyplot import imshow, axis, contour, contourf
import numpy as np
import os, sys

modeNum = 0
resultsDir = 'results'
normalizeSpatialGrowth = True   # Normalises spatial growth (to see absolute growth/decay)
ignoreTemporalGrowth = True     # Useful for visualisation

maxFactor = 0.5         # Plot contours from +/-maxFactor*max(abs(V))
P = 5                   # Number of time periods to plot movie over
Np = 20                 # Number of time steps per period
includeEndTime = False  # Do not include endtime if movie is going to loop
fps = 10

#vcodec = msmpeg4:vbitrate=800
vcodec = 'msmpeg4v2:vbitrate=800'
#vcodec = 'wmv2'

def frameZero():
    Np = 1
    P = 1
    #for fname in os.listdir(resultsDir):
    for fname in ['evals_R5000.mat']:
        plotTimes(fname)
    show()
    
def movie():
    #for fname in os.listdir(resultsDir):
    for fname in ['evals_R5000.mat']:
        plotTimes(fname,dumpFigs=True)
    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(int(fps)) + " -ovc lavc -lavcopts vcodec=" + vcodec + " -nosound -o " + fname + ".avi")
    os.system("rm -v *.png")

def plotTimes(fname,dumpFigs=False):

    vdict = loadmat('VARS.mat',struct_as_record=True)
    x = [vdict['pcx'][0][i] for i in xrange(len(vdict['pcx'][0]))]
    y = [vdict['pcy'][i][0] for i in xrange(len(vdict['pcy']))]
    x = array(x)
    y = array(y)
    L = float(vdict['LT'][0][0])

    # NOTE: figsize gives the figure size in inches @ 100dpi => 16=1200pixels
    fig = figure(figsize=(16,4))
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

        contourf(x,y,vmT,Z)
        axis('tight')
        
        if dumpFigs == True:
            fn = '_tmp%03d.png'%f
            print 'Saving frame', fn
            fig.savefig(fn)
        f += 1


