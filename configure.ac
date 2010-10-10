
AC_PREREQ(2.59)
AC_INIT([dolfin],[0.8.1-hpc])
AM_INIT_AUTOMAKE($PACKAGE_NAME, $PACKAGE_VERSION, no-define)
AM_MAINTAINER_MODE
AC_LANG(C++)

# Feature options.
AC_ARG_ENABLE(mpi,
              AC_HELP_STRING([--enable-mpi],
	      [Compile with support for MPI]),enable_mpi=yes, enable_mpi=no)

AC_ARG_ENABLE(openmp,
              AC_HELP_STRING([--enable-openmp],
	      [Compile with support for OpenMP]),enable_openmp=yes, enable_openmp=no)

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

if test $enable_mpi = yes; then
   AX_MPI	
   CXX="$MPICXX"
   LIBS="$MPILIBS $LIBS"
fi

if test $enable_openmp = yes; then
   AX_OPENMP
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
fi

# Checks for libraries.
AM_PATH_XML2( ,libxml2_found=yes, libxml2_found=no)
LIBS="$XML_LIBS $LIBS"

AX_CHECK_ZLIB
AX_BOOST_BASE([1.38],,AC_MSG_ERROR([suitable Boost not found]))

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])


#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile\
			  src/common/Makefile\
			  src/elements/Makefile\
			  src/fem/Makefile\
			  src/function/Makefile\
			  src/graph/Makefile\
			  src/io/Makefile])

AC_OUTPUT

 CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"