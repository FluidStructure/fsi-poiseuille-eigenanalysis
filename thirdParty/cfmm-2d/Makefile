# Macro Variables
# FMM files that mex functions depend upon
fmmfiles = treesolverlib kernellib calctreeupward calctreedownward treeparticlemerge boxtreebuild boxtree fmmboxclass containerlib
# Mex files that octave and matlab use
mexfiles = lambvortexfmm lambvortexvortfmm lambvortexnaive vortexelementnaive vortexelementfmm vortexsheetnaive
pyfiles = velocitymodule
# Resulting/needed files
fmmsources = $(addsuffix .cpp,${fmmfiles})
fmmobjects = $(addsuffix .o,${fmmfiles}) 
pyobj = $(addsuffix .o,${pyfiles})
octaveoutputs = $(addsuffix .mex,${mexfiles})
matlaboutputs = $(addsuffix .mexa64,${mexfiles})
pythonoutputs = $(addsuffix .so,${pyfiles})

# Module Install path
pythontoolbox = ../python_tools.git/

# Compilers used
CC = g++
LD = $(CC)
OCTAVEMEX = mkoctfile --mex
MATLABMEX = mex

# Complier and Linker flags
PYTHON = python2.6
LFLAGS = -shared -lm -lpthread -l${PYTHON}
CFLAGS = -pipe -O2 -fomit-frame-pointer -fpic -Wall -Wextra -D_REENTRANT -I /usr/include/${PYTHON}
#CFLAGS += -Weffc++
CFLAGS += -march=native

#------for intel core-2 chips (jarrad/mark/novak)
#CFLAGS += -msse4.1
#------for intel i7 chips (james/andrew/cluster-node/post-pro):
#CFLAGS += -msse4 -msse4.1 -msse4.2
#------for AMD systems
#CFLAGS += -msse3
#CFLAGS += -fno-inline -g

#------------------------------------------------------------------------

# List fake targets 
.PHONY: install uninstall python all_mex octave matlab clean usage
	
# Default rules...compile python and copy to ../pythontoolbox
# fix for broken matlab linking on multithreaded

install: python
	mkdir -p $(pythontoolbox)
	cp -f $(pythonoutputs) $(pythontoolbox)

uninstall:
	rm -rf $(pythontoolbox)*.so

debug: 
	@make clean
	@export DBF=-D_DEBUG_MODE;make -j4 -f Makefile install; export DBF=

python: $(pythonoutputs)
$(pythonoutputs):$(pyobj) $(fmmobjects)
	$(LD) $(LFLAGS) $< $(fmmobjects) -o $(pythonoutputs)

all_mex: octave
	make -j1 -f Makefile matlab

# Rest of rules	
octave: $(octaveoutputs)
%.mex:%.cpp $(fmmobjects)
	@echo "WARNING: Octave interface still works, but support ceased on 08-2009, python is the only active port...Press <enter> to continue"
	@read tmp
	$(OCTAVEMEX) $< $(fmmobjects)

matlab: $(matlaboutputs)	
%.mexa64:%.cpp $(fmmobjects)
	@echo "WARNING: Matlab support has been deprecated since 08-2009...Press <enter> to continue"
	@read tmp
	$(MATLABMEX) $< $(fmmobjects)

%.o:%.cpp 
	$(CC) $(CFLAGS) $(DBF) -c $<

fmmobjects: $(fmmobjects)

testfast:
	rm -rf *.o
	$(CC) $(CFLAGS) -c testfmm.cpp testfmmmerge.cpp $(fmmsources)
	$(LD) -lm -lpthread testfmm.o $(fmmobjects) -o testfmm.bin
	$(LD) -lm -lpthread testfmmmerge.o $(fmmobjects) -o testfmmmerge.bin

test:
	rm -rf *.o
	$(CC) $(CFLAGS) -fno-inline -g -c testfmm.cpp testfmmmerge.cpp $(fmmsources)
	$(LD) -lm -lpthread testfmm.o $(fmmobjects) -o testfmm.bin
	$(LD) -lm -lpthread testfmmmerge.o $(fmmobjects) -o testfmmmerge.bin


usage: 
	@echo "'make usage' display this message"
	@echo ""
	@echo "'make install' to build and install python fmm module"
	@echo ""
	@echo "'make uninstall' to remove installed python module"
	@echo ""
	@echo "'make python' to build python module only"
	@echo ""
	@echo "'make all_mex' to build octave and matlab files"
	@echo ""
	@echo "'make python' to build python files only"
	@echo ""
	@echo "'make octave' to build octave files only"
	@echo ""
	@echo "'make fmmobjects' to build fmm objects only"
	@echo ""
	@echo "'make matlab' to build matlab files only"
	@echo ""
	@echo "'make clean' to remove all built files"
clean:
	rm -rf *.o *.me* a.out *.so
