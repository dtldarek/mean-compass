CC=clang
CXX=clang++
RM=rm -f
CPPFLAGS=-g -std=c++11 -I/usr/include/eigen3/
LDFLAGS=-g
LDLIBS=-lmpfr -lgmp

SRCS=mean-compass.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: mean-compass

mean-compass: $(OBJS)
	$(CXX) $(LDFLAGS) -o mean-compass $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) *~ .depend

include .depend
