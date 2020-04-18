# VERSION    = core or simp
# SOLVERNAME = name of the SAT solver
# SOLVERDIR  = subdirectory of the SAT solver
# NSPACE     = namespace of the SAT solver
#
# e.g. minisat compilation with core version:
#
# VERSION    = core
# SOLVERNAME = "Minisat"
# SOLVERDIR  = minisat
# NSPACE     = Minisat
#
VERSION    = core
SOLVERNAME = "Glucose4.1"
SUPERSOLVERNAME=TT-Open-WBO-Inc#TT-Open-WBO-Inc Loandra Open-WBO-Inc LinSBPS SATLike
SUPERSOLVERNAMEID=1#1 2 3 4 5

NSPACE     = Glucose
ifeq ($(SUPERSOLVERNAMEID), 5)
Dist: solver/SATLike/basis_pms.h solver/SATLike/pms.h solver/SATLike/pms.cpp rapidjson/*.h rapidjson/msinttypes/*.h rapidjson/internal/*.h rapidjson/error/*.h problem/*.h
	g++ -std=c++11 main.cc -DMAXSATNID=$(SUPERSOLVERNAMEID)  -O3  -o timetabler
endif
ifneq ($(SUPERSOLVERNAMEID), 5)
SOLVERDIR  = solver/$(SUPERSOLVERNAME)/solvers/glucose4.1
# THE REMAINING OF THE MAKEFILE SHOULD BE LEFT UNCHANGED
EXEC       = timetabler
DEPDIR     = mtl utils core
DEPDIR     +=  ../../../$(SUPERSOLVERNAME) ../../encodings ../../algorithms ../../graph ../../classifier ../../clusterings ../../../../problem   ../../../../rapidXMLParser
MROOT      = $(PWD)/$(SOLVERDIR)
LFLAGS     += -lgmpxx -lgmp
CFLAGS     =  -DMAXSATNID=$(SUPERSOLVERNAMEID)  -O3 -Wall -Wno-parentheses -std=c++11 -DNSPACE=$(NSPACE) -DSOLVERNAME=$(SOLVERNAME) -DVERSION=$(VERSION)
ifeq ($(VERSION),simp)
DEPDIR     += simp
CFLAGS     += -DSIMP=1
ifeq ($(SOLVERDIR),glucored)
LFLAGS     += -pthread
CFLAGS     += -DGLUCORED
DEPDIR     += reducer glucored
endif
endif
include $(MROOT)/mtl/template.mk
endif
