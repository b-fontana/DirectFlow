CC = g++
EXEC = v1_beam.exe
RM = rm -f

DEBUG_LEVEL     = -g -fdiagnostics-color=always
EXTRA_CCFLAGS   = -Wall -std=c++17 -O -pedantic -pedantic-errors -Wformat -Wformat=2 \
	-Wformat-nonliteral -Wformat-security  \
        -Wformat-y2k \
        -Wimport  -Winit-self  -Winline \
        -Winvalid-pch   \
        -Wunsafe-loop-optimizations -Wmissing-braces \
        -Wmissing-field-initializers -Wmissing-format-attribute   \
        -Wmissing-include-dirs -Wmissing-noreturn
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS         = $(CXXFLAGS)
ROOTFLAGS = `root-config --cflags --ldflags --evelibs` -lMinuit -lrt
OPENGLFLAGS = -L/usr/lib/x86_64-linux-gnu/ -lGL -lGLX -lGLdispatch

BASEDIR := $(shell pwd)
SRCDIR  := src
INCDIR  := include

SRCS := $(wildcard *.cc) \
	$(wildcard $(SRCDIR)/*.cc)

OBJS := $(patsubst %.cc, %.o, $(SRCS))
## Remove "special" sources
# SPECIAL_SOURCES := divisor.cc multiplier.cc
## SOURCES := $(filter-out $(SPECIAL_SOURCES),$(SOURCES))

.PHONY: all clean
.DEFAULT_GOAL = all

all: $(EXEC)

$(EXEC): $(OBJS) ;
	$(CC) $(CCFLAGS) $^ $(ROOTFLAGS) $(OPENGLFLAGS) -o $@
	@echo Executable $(EXEC) created.

%.o: %.cc
	$(CC) $(CCFLAGS) -c $< $(ROOTFLAGS) -I$(BASEDIR) -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.cc
	$(CC) $(CCFLAGS) -c $< $(ROOTFLAGS) -I$(BASEDIR) -o $@

clean:
	$(RM) $(OBJS) $(EXEC)
