check_PROGRAMS = io_test
io_test_SOURCES = io_test.cpp
exe_CPPFLAGS= -I@top_srcdir@/include @CHECK_CFLAGS@ $(CPPFLAGS)
exe_LDADD = $(top_builddir)/src/libdolfin.la
exe_LDFLAGS = $(LDFLAGS) $(LIBS) @CHECK_LIBS@
