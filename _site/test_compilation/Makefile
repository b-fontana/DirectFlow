CC = g++
EXEC = v1_beam.exe
RM = rm -r

BASEDIR := $(shell pwd)
SRCDIR  := src
INCDIR  := include
DEPDIR  := .deps

DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
DEBUG_LEVEL := -g -fdiagnostics-color=always
EXTRA_CCFLAGS := -Wall -std=c++17 -O -pedantic -pedantic-errors
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS         = $(CXXFLAGS)

SRCS := $(wildcard *.cc) \
	$(wildcard $(SRCDIR)/*.cc)

OBJS := $(patsubst %.cc, %.o, $(SRCS))

DEPFILES := $(patsubst %.cc, $(DEPDIR)/%.d, $(notdir $(SRCS)))

.PHONY: all clean
.DEFAULT_GOAL = all

all: $(DEPDIR) $(EXEC)

$(EXEC): $(OBJS)
	@echo 3
	@echo Dependencies: $(wildcard $(DEPFILES))
	$(CC) $(CCFLAGS) $^ -o $@
	@echo Executable $(EXEC) created.

%.o: %.cc
%.o: %.cc Makefile
	@echo 1	
	$(CC) $(DEPFLAGS) $(CCFLAGS) -c $< -I$(BASEDIR) -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.cc $(DEPDIR)/%.d | $(DEPDIR)
	@echo 2
	$(CC) $(DEPFLAGS) $(CCFLAGS) -c $< -I$(BASEDIR) -o $@

$(DEPDIR):
	mkdir -p $@

$(DEPFILES):

clean:
	$(RM) $(OBJS) $(EXEC) $(DEPDIR)

-include $(wildcard $(DEPFILES))
