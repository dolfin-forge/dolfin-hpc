SUBDIRS=petsc

AM_CPPFLAGS = -I@top_srcdir@/include

noinst_LTLIBRARIES = libdolfin_la.la

libdolfin_la_la_LIBADD = petsc/libdolfin_la_petsc.la

libdolfin_la_la_SOURCES = DefaultFactory.cpp\
			                    JANPACKFactory.cpp\
			                    JANPACKKrylovSolver.cpp\
			                    JANPACKAMGSolver.cpp\
			                    JANPACKMat.cpp\
			                    JANPACKVec.cpp\
			                    KrylovSolver.cpp\
			                    LinearSolver.cpp\
			                    LAPACKMatrix.cpp\
			                    LAPACKVector.cpp\
			                    SparsityPattern.cpp\
			                    TrilinosVector.cpp
