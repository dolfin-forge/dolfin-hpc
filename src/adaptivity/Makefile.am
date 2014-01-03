AM_CPPFLAGS = -I@top_srcdir@/include
noinst_LTLIBRARIES = libdolfin_adaptivity.la
libdolfin_adaptivity_la_SOURCES = AdaptiveRefinement.cpp\
				SpaceTimeFunction.cpp

if ENABLE_UFL
UFCDIR= ufc2
else
UFCDIR= ufc1
endif

AM_CPPFLAGS+= -I@top_srcdir@/include/dolfin/adaptivity/${UFCDIR}
libdolfin_adaptivity_la_SOURCES+= ${UFCDIR}/AdaptiveRefinementProjectScalar.cpp\
				${UFCDIR}/AdaptiveRefinementProjectVector.cpp
