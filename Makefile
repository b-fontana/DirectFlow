CC = g++
EXEC = v1_beam.exe
RM = rm

DEBUG_LEVEL     = -g
EXTRA_CCFLAGS   = -Wall
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS         = $(CXXFLAGS)
ROOTFLAGS = `root-config --cflags --ldflags --evelibs` -lMinuit -lrt \
	-I/home/bruno/miniconda3/envs/DirectFlow/include/
OPENGLFLAGS = -L/usr/lib/x86_64-linux-gnu/ -lGL -lGLX -lGLdispatch

SRCS  = $(wildcard *.cc)
OBJS := $(patsubst %.cc, %.o, $(SRCS))
## Remove "special" sources
# SPECIAL_SOURCES := divisor.cc multiplier.cc
## SOURCES := $(filter-out $(SPECIAL_SOURCES),$(SOURCES))

.PHONY: all clean
.DEFAULT_GOAL = all

all: ${EXEC}

%.o: %.cc
	$(CC) -c $(CCFLAGS) $< $(ROOTFLAGS)

$(EXEC): $(OBJS) ;
	$(CC) $(CCFLAGS) $^ $(ROOTFLAGS) $(OPENGLFLAGS) -o $@

clean:
	$(RM) $(OBJS) $(EXEC)
