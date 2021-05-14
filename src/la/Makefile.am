SUBDIRS=janpack petsc trilinos

AM_CPPFLAGS = -I@top_srcdir@/include

noinst_LTLIBRARIES = libdolfin_la.la

libdolfin_la_la_LIBADD = janpack/libdolfin_la_janpack.la\
                         petsc/libdolfin_la_petsc.la\
                         trilinos/libdolfin_la_trilinos.la

libdolfin_la_la_SOURCES = DefaultFactory.cpp\
                          KrylovSolver.cpp\
                          LinearSolver.cpp\
                          LAPACKMatrix.cpp\
                          LAPACKVector.cpp\
                          SparsityPattern.cpp
