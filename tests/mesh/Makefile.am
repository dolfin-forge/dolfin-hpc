check_PROGRAMS = exe
exe_SOURCES = test.cpp
exe_CPPFLAGS = -I@top_srcdir@/include @CHECK_CFLAGS@ $(CPPFLAGS)
exe_LDADD = $(top_builddir)/src/libdolfin.la
exe_LDFLAGS = $(LDFLAGS) $(LIBS) @CHECK_LIBS@


EXTRA_DIST = CellTypes/CellTypes.h\
             SubDomain/SubDomain.h\
             VertexNormal/VertexNormal.h\
             StructuredGrid/StructuredGrid.h\
             UnitMeshes/UnitMeshes.h

