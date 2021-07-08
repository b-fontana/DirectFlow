DEBUG_LEVEL     = -g
EXTRA_CCFLAGS   = -Wall
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS         = $(CXXFLAGS)

ROOTFLAGS = `root-config --cflags --ldflags --evelibs` -lMinuit -lrt -I/home/bruno/miniconda3/envs/DirectFlow/include/
OPENGLFLAGS = -L/usr/lib/x86_64-linux-gnu/ -lGL -lGLX -lGLdispatch

CC = g++

EXEC = v1_beam.exe

RM = rm
# g++ -g v1_beam.cc tracking.cc `root-config --cflags --ldflags --evelibs` -lMinuit -lrt -I/home/bruno/miniconda3/envs/DirectFlow/include/ -L/usr/lib/x86_64-linux-gnu/ -lGL -lGLX -lGLdispatch -o "$(EXE)";
# echo "Executable $(EXEC) ready."

SRCS  = $(wildcard *.cc)
## Remove "special" sources
# SPECIAL_SOURCES := divisor.cc multiplier.cc
## SOURCES := $(filter-out $(SPECIAL_SOURCES),$(SOURCES))
OBJS := $(patsubst %.cc, %.o, $(SRCS))

.PHONY: all clean
.DEFAULT_GOAL = all

all: ${EXEC} ;

%.o: %.cc ;
	gcc -c $(CCFLAGS) $<;

$(EXEC): $(OBJS) ;
	$(CC) $(CCFLAGS) $^ $(ROOTFLAGS) $(OPENGLFLAGS) -o $@;

clean:
	$(RM) $(OBJS) $(EXEC);
