# Define C++ compiler
CXX = g++

# Define linker
LINKER = g++

# Debug flag: either nothing or -ggdb
DEBUGFLAG = 

# Opitmization flag
OPTFLAG = -O2

# Path to StPicoEvent (either absolute or relative)
FEMTOPATH = /home/kirill/Soft/StFemtoEvent
REFMULTPATH = /home/kirill/Soft/StRefMultCorr

# ROOT flag: -D__ROOT__
ROOTFLAG =

# Define C++ flags (-D_VANILLA_ROOT_ is needed to avoid StMessMgr confusion)
# -std=c++11
CXXFLAGS = $(shell root-config --cflags) -I. -I$(FEMTOPATH) -I$(REFMULTPATH) $(OPTFLAG) -fPIC -Wall -pipe
CXXFLAGS += -D_VANILLA_ROOT_ $(DEBUGFLAG)

# Define libraries to be added
LIBS = $(shell root-config --libs) 
LIBS += -L$(FEMTOPATH) -lStFemtoDst
LIBS += -L$(REFMULTPATH) -lStRefMultCorr

# Define path to include directories
INCS = $(shell root-config --incdir)

# Linker flags //  -D__ROOT__
LDFLAGS = $(shell root-config --cflags) -I. -I$(FEMTOPATH) -I$(REFMULTPATH) $(OPTFLAG) -fPIC -Wall -pipe -D_VANILLA_ROOT_

all: femtoAnalyzerStandalone

# Compile C++ files
femtoAnalyzerStandalone: femtoAnalyzerStandalone.o
	$(LINKER) $< $(LDFLAGS) $(LIBS) $(GLIBS) -o $@

# Compile C++ files
%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo 'Compiling: $@ $<'

# PHONY options
.PHONY: clean distclean

# Remove object files and dictionary
clean:
	rm -vf *.o

# Remove object files, dictionary and library
distclean:
	rm -vf *.o femtoAnalyzerStandalone
