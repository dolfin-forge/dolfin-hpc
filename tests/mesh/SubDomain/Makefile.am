TESTS = exe
check_PROGRAMS = exe
exe_SOURCES = test.cpp
exe_CPPFLAGS= -I@top_srcdir@/include $(CPPFLAGS)
exe_LDFLAGS = $(top_builddir)/src/.libs/libdolfin.a $(LDFLAGS) $(LIBS) @CHECK_LIBS@
