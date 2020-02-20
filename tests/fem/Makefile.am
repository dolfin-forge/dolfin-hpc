check_PROGRAMS = exe
exe_SOURCES = test.cpp
nodist_exe_SOURCES = ../tests.h
exe_CPPFLAGS= -I@top_builddir@/tests -I@top_srcdir@/include @CHECK_CFLAGS@ $(CPPFLAGS)
exe_LDADD = $(top_builddir)/src/libdolfin.la
exe_LDFLAGS = $(LDFLAGS) $(LIBS) @CHECK_LIBS@


EXTRA_DIST = BoundaryNormal/BoundaryNormal.h\
      	     DiscreteSpaces/DiscreteSpaces.h\
            DofMap/DofMap.h\
            DofNumbering/DofNumbering.h\
      	     FiniteElement/FiniteElement.h\
      	     FiniteElementSpace/FiniteElementSpace.h\
            Form/Form.h
