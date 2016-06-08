RM=rm -f
CXXFLAGS=-Wall -Wextra -Wpedantic \
		 -Wshadow -Wpointer-arith -Wcast-qual \
		 -Wstrict-overflow=2 -Wlogical-op -Wfloat-equal \
		 -Wswitch-enum -Wswitch-default -Wredundant-decls \
		 -fstrict-overflow \
		 -MMD -std=c++11 -isystem/usr/include/eigen3/
LDFLAGS=-Wl,-z,now
LDLIBS=-lmpfr -lgmp -lboost_random -lboost_program_options

SUBDIRS=.depends build
SRCS=src/utils.cc src/utf8_io.cc src/mean_compass.cc
OBJS=$(SRCS:src/%.cc=build/%.o)

all: CXXFLAGS += -DNDEBUG -O3
all: mean-compass

perf: CXXFLAGS += -DNDEBUG -O3 -g \
	              -fno-omit-frame-pointer \
				  -fno-inline-functions \
				  -fno-inline-functions-called-once \
				  -fno-optimize-sibling-calls
perf: mean-compass

debug: CXXFLAGS += -DDEBUG -g -Og
debug: LDFLAGS += -g
debug: mean-compass

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
