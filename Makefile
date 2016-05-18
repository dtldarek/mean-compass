CC=gcc
CXX=g++
RM=rm -f
CXXFLAGS=-g -Wall -MMD -std=c++11 -I/usr/include/eigen3/
LDFLAGS=-g
LDLIBS=-lmpfr -lgmp -lboost_random -lboost_program_options

SUBDIRS=.depends build
SRCS=src/mean_compass.cc src/utils.cc
OBJS=$(SRCS:src/%.cc=build/%.o)

all: mean-compass

mean-compass: $(OBJS)
	@echo Linking $@
	$(CXX) $(LDFLAGS) -o mean-compass $(OBJS) $(LDLIBS) 

build/%.o: src/%.cc
	@echo Compiling $@
	@mkdir -p $(SUBDIRS)
	$(CXX) $(CXXFLAGS) -MMD -MF $(patsubst build/%.o,.depends/%.d,$@) -c -o $@ $<

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) *~ .depends/*.d

-include $(SRCS:src/%.cc=.depends/%.d)
