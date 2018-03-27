check_PROGRAMS = exe
exe_SOURCES = test.cpp
nodist_exe_SOURCES = ../tests.h
exe_CPPFLAGS= -I@top_builddir@/tests -I@top_srcdir@/include @CHECK_CFLAGS@ $(CPPFLAGS)
exe_LDADD = $(top_builddir)/src/libdolfin.la
exe_LDFLAGS = $(LDFLAGS) $(LIBS) @CHECK_LIBS@


EXTRA_DIST = Cell/UFL_Cell.h\
             Domain/UFL_Domain.h\
             FiniteElement/UFL_FiniteElement.h\
      	     MixedElement/UFL_MixedElement.h\
             Space/UFL_Space.h\
      	     type/UFL_type.h