##  Template makefile for Standard, Debug, Release, and Release-static versions
##
##    "make"    for the standard version (optimized, but with debug information and assertions active)
##    "make r"  for a release version.
##    "make rs" for a statically linked release version.
##    "make d"  for a debug version (no optimizations).
##
##  Pace-specific versions (same as release with custom flags)
##    "make heuristic"
##    "make exact"
##    "make cutwidth"
##
##    "make lite" (a very fast version)

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
HCOBJS     = $(addsuffix h,  $(COBJS))
ECOBJS     = $(addsuffix e,  $(COBJS))
CCOBJS     = $(addsuffix c,  $(COBJS))
LCOBJS     = $(addsuffix l,  $(COBJS))

CXX       ?= clang++

CFLAGS    ?= -Wall -Wno-deprecated-declarations -std=c++17
LFLAGS    ?= -Wall -flto -fuse-ld=lld

COPTIMIZE ?= -O3

CFLAGS    += -I$(MROOT)

.PHONY : s d r rs heuristic exact cutwidth lite clean

## Executable names
s:	       $(EXEC)
d:	       $(EXEC)_debug
r:	       $(EXEC)_release
rs:	       $(EXEC)_static
heuristic: $(EXEC)_heuristic
exact:	   $(EXEC)_exact
cutwidth:	 $(EXEC)_cutwidth
lite:	     $(EXEC)_lite

## Compile options
%.o:	CFLAGS +=$(COPTIMIZE) -g -D DEBUG
%.od:	CFLAGS +=-O0 -g -D DEBUG
%.or:	CFLAGS +=$(COPTIMIZE) -g -D NDEBUG
%.oh:	CFLAGS +=$(COPTIMIZE) -g -D NDEBUG -D HEURISTIC
%.oe:	CFLAGS +=$(COPTIMIZE) -g -D NDEBUG -D EXACT
%.oc:	CFLAGS +=$(COPTIMIZE) -g -D NDEBUG -D CUTWIDTH
%.ol:	CFLAGS +=$(COPTIMIZE) -g -D NDEBUG -D LITE

## Link options
$(EXEC):		       LFLAGS += -g
$(EXEC)_debug:		 LFLAGS += -g
$(EXEC)_static:		 LFLAGS += --static
$(EXEC)_release:	 LFLAGS +=
$(EXEC)_heuristic: LFLAGS +=
$(EXEC)_exact:	   LFLAGS +=
$(EXEC)_cutwidth:	 LFLAGS +=
$(EXEC)_lite:	     LFLAGS +=

## Dependencies
$(EXEC):		       $(COBJS)
$(EXEC)_debug:		 $(DCOBJS)
$(EXEC)_release:	 $(RCOBJS)
$(EXEC)_static:		 $(RCOBJS)
$(EXEC)_heuristic: $(HCOBJS)
$(EXEC)_exact:     $(ECOBJS)
$(EXEC)_cutwidth:  $(CCOBJS)
$(EXEC)_lite:  	   $(LCOBJS)

## Build rule
%.o %.od %.or %.oh %.oe %.oc %.ol:	%.cpp
	@echo Compiling: $(subst $(MROOT)/,,$@)
	@$(CXX) $(CFLAGS) -c -o $@ $<

## Linking rules (standard/profile/debug/release)
$(EXEC) $(EXEC)_debug $(EXEC)_release $(EXEC)_static $(EXEC)_heuristic $(EXEC)_exact $(EXEC)_cutwidth $(EXEC)_lite:
	@echo Linking: "$@ ( $(foreach f,$^,$(subst $(MROOT)/,,$f)) )"
	@$(CXX) $^ $(LFLAGS) -o $@

## Clean rule
clean:
	@rm -f $(EXEC) $(EXEC)_* $(COBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) $(HCOBJS) $(ECOBJS) $(CCOBJS) $(LCOBJS) \
	       log_* *.dimacs.gz depend.mk tmp*.res

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
