AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(PH_HOME)/install/include

lib_LTLIBRARIES = \
  libanalyzePhotons.la 

noinst_HEADERS = \
  analyzePhotonsLinkDef.h 

include_HEADERS = \
  analyzePhotons.h \
  plotPhotons.h \
  deadMap.h \
  dchDeadMap.h \
  emcalDeadMap.h

libanalyzePhotons_la_SOURCES = \
  analyzePhotons.cpp \
  plotPhotons.cpp \
  deadMap.cpp \
  dchDeadMap.cpp \
  emcalDeadMap.cpp \
  analyzePhotons_Dict.C

libanalyzePhotons_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(LD_LIBRARY_PATH) \
  -L$(ROOTSYS)/lib/libRIO

# Rule for generating CINT dictionaries from class headers.
analyzePhotons_Dict.C: \
  analyzePhotons.h \
  plotPhotons.h \
  deadMap.h \
  dchDeadMap.h \
  emcalDeadMap.h \
  analyzePhotonsLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
