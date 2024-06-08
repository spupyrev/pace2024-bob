##  Template makefile for Standard, Debug, Release, and Release-static versions
##
##    "make"    for the standard version (optimized, but with debug information and assertions active)
##    "make r"  for a release version.
##    "make rs" for a statically linked release version.
##    "make d"  for a debug version (no optimizations).

EXEC      = pace
DEPDIR    = src
MROOT = $(PWD)

PWD        = $(shell pwd)
EXEC      ?= $(notdir $(PWD))

CSRCS      = $(wildcard $(PWD)/*.cpp)
DSRCS      = $(foreach dir, $(DEPDIR), $(filter-out $(MROOT)/$(dir)/main*.cpp, $(wildcard $(MROOT)/$(dir)/*.cpp)))
CHDRS      = $(wildcard $(PWD)/*.h)
COBJS      = $(CSRCS:.cpp=.o) $(DSRCS:.cpp=.o)

PCOBJS     = $(addsuffix p,  $(COBJS))
DCOBJS     = $(addsuffix d,  $(COBJS))
RCOBJS     = $(addsuffix r,  $(COBJS))

CXX       ?= clang++
# CXX       ?= g++

CFLAGS    ?= -Wall -Wno-deprecated-declarations -std=c++17
LFLAGS    ?= -Wall -flto -fuse-ld=lld # remove lld option for g++

COPTIMIZE ?= -O3

CFLAGS    += -I$(MROOT) -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS

.PHONY : s d r rs clean

s:	$(EXEC)
d:	$(EXEC)_debug
r:	$(EXEC)_release
rs:	$(EXEC)_static

## Compile options
%.o:			CFLAGS +=$(COPTIMIZE) -g -D DEBUG
%.od:			CFLAGS +=-O0 -g -D DEBUG
%.or:			CFLAGS +=$(COPTIMIZE) -g -D NDEBUG

## Link options
$(EXEC):		      LFLAGS += -g
$(EXEC)_debug:		LFLAGS += -g
$(EXEC)_release:	LFLAGS +=
$(EXEC)_static:		LFLAGS += --static

## Dependencies
$(EXEC):		$(COBJS)
$(EXEC)_debug:		$(DCOBJS)
$(EXEC)_release:	$(RCOBJS)
$(EXEC)_static:		$(RCOBJS)

## Build rule
%.o %.od %.or:	%.cpp
	@echo Compiling: $(subst $(MROOT)/,,$@)
	@$(CXX) $(CFLAGS) -c -o $@ $<

## Linking rules (standard/profile/debug/release)
$(EXEC) $(EXEC)_debug $(EXEC)_release $(EXEC)_static:
	@echo Linking: "$@ ( $(foreach f,$^,$(subst $(MROOT)/,,$f)) )"
	@$(CXX) $^ $(LFLAGS) -o $@

## Clean rule
clean:
	@rm -f $(EXEC) $(EXEC)_debug $(EXEC)_release $(EXEC)_static \
	  $(COBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) log_* *.dimacs.gz depend.mk tmp.res

## Make dependencies
depend.mk: $(CSRCS) $(CHDRS)
	@echo Making dependencies
	@$(CXX) $(CFLAGS) -I$(MROOT) \
	   $(CSRCS) -MM | sed 's|\(.*\):|$(PWD)/\1 $(PWD)/\1r $(PWD)/\1d $(PWD)/\1p:|' > depend.mk
	@for dir in $(DEPDIR); do \
	      if [ -r $(MROOT)/$${dir}/depend.mk ]; then \
		  echo Depends on: $${dir}; \
		  cat $(MROOT)/$${dir}/depend.mk >> depend.mk; \
	      fi; \
	  done
