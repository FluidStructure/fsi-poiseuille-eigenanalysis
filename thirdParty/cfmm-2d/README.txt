//==============|    
//  NAME        : README.txt
//  AUTHOR      : Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
//  MODIFIED    : 16.06.2009
//  DESCRIPTION : information relating to installation and other stuff about the c-FMM. 
//==============|

-------------------------------------------------------------------
INSTALLATION-INSTRUCTIONS:
At the shell "cd" to the fmm directory then type:
    make usage; # to show make options
To try C benchmark case do:
    make testfast;
To setup python bindings do:
    make clean; make -j5 python;
If using with the python_toolbox:
    make clean; make; # equivalent to "make python; make install"


Then all should be compiled automatically with some fairly aggressive compiler optimizations set in MAKEFILE for an intel core 2 duo system. 
Remove the -03 and -march=core2 for safer build flags, possibly replace -fpic with -fPIC also.
-------------------------------------------------------------------


-------------------------------------------------------------------
USE FORMAT:
    [u,v] = lambVortexNaive(precision, threads, x,y,str,coresqrd,ex,ey); 
    [u,v] = lambvortexfmm( maxPvalue, targetsPerBox, threads, x,y,str,coresqrd,ex,ey);    
Where u and v returned are the velocities at each target inputted at ex,ey with particles/source info x,y,str,core.
    
    [u,v] = vortexelementfmm(xl,yl,strl, xr, yr, strr, evalx, evaly, accuracy, targetsPerBox, threads);
    [u,v] = vortexelementnaive(xl,yl,strl, xr, yr, strr, evalx, evaly, accuracy, threads);
Where u and v returned are the velocities at each target inputted at ex,ey with each sheet being described by its too ends xl,yl,strl, and xr, yr, strr

    [u,v] = vortexsheetnaive(xl,yl,xr, yr, strr, evalx, evaly, threads);
Where u and v returned are the velocities at each target inputted at ex,ey with each sheet being described by its too ends xl,yl,and xr, yr, strr


USE EXAMPLE:
In octave to use the fmm function, to try a N^2 benchmark for a dualcore try;
    maxPvalue=9; targetsPerBox=15; n=1000000; coreSizeSqrd=0; x=randn(n,1); y=randn(n,1); str=randn(n,1); core=str*0+coreSizeSqrd; ex=x; ey=y; for threads=4:4; for run=1:1; disp(threads);  tic; [u,v] = lambvortexfmm( maxPvalue, targetsPerBox, threads, x,y,str,core,ex,ey); time(threads,run)=toc; disp(time(threads,run));end;end;
    
    threads=4; maxPvalue=9; targetsPerBox=5; n=500000; x=randn(n,1); y=randn(n,1); str=randn(n,1); str2=randn(n,1); len=0.00025/2; ex=x; ey=y; tic; [u,v] = vortexelementfmm( x-len, y, str, x+len, y, str, ex, ey, maxPvalue, targetsPerBox, threads); disp(toc);

In python (from within cfmm install directory only):
    python -c "import random, time, velocitymodule; maxPvalue=9; targetsPerBox=15; n=1000000; coreSizeSqrd=0; x=[random.random() for i in xrange(n)]; y=[random.random() for i in xrange(n)]; s=[random.random() for i in xrange(n)]; core=[coreSizeSqrd]*n; ex=x; ey=y; threads=4; run=1; print threads; tic=time.time(); ou =[0.0]*n; ov=[0.0]*n; velocitymodule.lamb_vortex(x,y,s,core, ex, ey, ou,ov, threads, maxPvalue, targetsPerBox); print 'time taken' , str(time.time() - tic)"

In c:
    make testfast; ./a.out

If python_tools and cfmm are installed together, from within python_tools directory do:
    python -c "import random, time, velocitylib; n=1000000; coreSizeSqrd=0; x=[random.random() for i in xrange(n)]; y=[random.random() for i in xrange(n)]; s=[random.random() for i in xrange(n)]; core=[coreSizeSqrd]*n; ex=x; ey=y; print 'Started with Defaults; threads=',4, 'fmm_accuracy=', 0.0001, 'fmm_targetsperbox=', 15  ; tic=time.time(); (ou,ov) = velocitylib.lamb_vortex_vel(x,y,s,core, ex, ey); print 'time taken' , str(time.time() - tic)"

-------------------------------------------------------------------
Current performance (22.10.08) coreSizeSqrd=0 p=9, ppb=15, threads=4:
    ~0.5s sim: 
        1.8ghz C2D laptop = 15,000 particles
        3.0ghz C2Q UniPc  = 40,000 particles
        
    1,000,000 particles:
        1.8ghz C2D laptop = 34s
        3.7ghz C2Q UniPc  = 8.5s
        
-------------------------------------------------------------------   

KNOWN BUGS: 
    Matlab whinges "version `GLIBCXX_3.4.9' not found"
    * for fix see " http://brahms.pbwiki.com/Errors " , second from bottom
    -   on jarrads 64bit system it invloves relocating the 3 offending files from /home/matlab2008a/sys/os/glnxa64/ to a new home
   
   cmds: 
   cd /home/matlab20*
   cd sys/os/glnxa64/ ; mkdir old; mv -v libgcc_s.so.1 libstdc++.so.6 libstdc++.so.6.0.8 old/ 

    WARNING: MATLAB IMPLEMENTATION DEPRECATED AND TAKES IN CORESQRD AND STR/2PI
-------------------------------------------------------------------        

EXTENSION REQUIREMENTS AND PROCESS
Should the need arrise to implement a new function for the purpose of use in matlab, the required modifications are (in a somewhat backwards order);
    - create a matlab/octave .m file to interface into the generated *.mex file into matlab 
    - edit Makefile to build your new c function + corresponding mex files
    - create the gateway mex function (in C) that matlab will use to interface with the new treesolverlib module
    - declare the new treesolverlib module in its header file treesolverlib.hpp and then implement it into the treesolverlib.cpp file.
        - implement the required virtual/abstract functions to complete module as:
        - transformFMMvel how to convert from the original fmm algorithm to desired particle in the far field.
        - kernelVel function...if desired implement into kernellib (.hpp + .cpp)
    Examles types:  farfield(fmm) + nearfield   -> lambvortexfmm
                    nearfield only              -> lambvortexvortfmm
                    farfield(direct) + nearfield-> lambvortexnaive

-------------------------------------------------------------------

FUTURE IMPROVEMENTS/FEATURES (Assigned: Jarrad):
    *   Implement GPU programming (MPI functionality not required).

VERSION CHAGELOG:
    Version 1.1.0 (20.04.10)
        *   Fixed memory leak in python interface files.
        *   Included lamb and panel naive* codes. (naive in accuracy, not implementation)
        *   Most function pre-multiply by 2Pi and coresize.
        *   Lambvortexvort function added.
        *   Added lambvortexmerge FMM function.

    Version 1.0.9 (02.10.09)
        *   Improved memory efficiency and speed.
        *   Included new test case.

    Version 1.0.8 (16.06.09)  
        *   Implemented lambVortexVortFmm to calculate vorticity fields.

    Version 1.0.7 (16.03.09)  
        *   Simplified fmm boxes.
        *   Included libhoard possibility.
        *   Changed lamb vortex kernels to take in sigma squared.
        
    Version 1.0.6 (20.10.08)  
        *   Implemented naive method for semi-infinite sheets.

    Version 1.0.5 (07.10.08)  
        *   Fully multi-threaded upward and buildtree steps.
        *   Implemented linear/constant discrete vortex sheets.
        
    Version 1.0.4 (04.09.08)    
        *   Implemented polymorphism and virtual functions.
        *   Included function for a multithreaded naive vortex code.
        *   Included switches that allow the far field effects to be turned off, including allowing the tree to still work at level 0.
        *   Implemented targets and sources seperately.
    
    Version 1.0.3 (28.08.08)
        *   Implemented dynamic-R-size optimization to the FMM algorithm.
        *   Removed extra safety checks and segfault proofing.
        *   Inlined many simple functions.

    Version 1.0.2 (26.08.08)
        *   Included posix multithreaded (pthread) support.

    Version 1.0.1 (25.08.08)
        *   Fixed major memory leaks.
        *   Fixed problem that it wouldnt use inputted PPB + maxtree values.
        *   Support added for lamb vortices.
