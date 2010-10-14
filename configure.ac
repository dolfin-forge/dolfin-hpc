
AC_PREREQ(2.59)
AC_INIT([dolfin],[0.8.1-hpc])
AC_DEFINE_UNQUOTED(DOLFIN_VERSION,["${PACKAGE_VERSION}"],[Version of DOLFIN])
AM_INIT_AUTOMAKE($PACKAGE_NAME, $PACKAGE_VERSION, no-define)
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
# Feature options.
AC_ARG_ENABLE(mpi,
              AC_HELP_STRING([--enable-mpi],
	      [Compile with support for MPI]),[enable_mpi=${enableval}], [enable_mpi=no])
AC_ARG_ENABLE(openmp,
              AC_HELP_STRING([--enable-openmp],
	      [Compile with support for OpenMP]),[enable_openmp=${enableval}], [enable_openmp=no])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET

if test "x${enable_mpi}" = xyes; then
   AX_MPI	
   CXX="$MPICXX"
   LIBS="$MPILIBS $LIBS"
fi

if test "x${enable_openmp}" = xyes; then
   AX_OPENMP
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
fi

# Checks for libraries.
AM_PATH_XML2(,,AC_MSG_ERROR([libxml2 not found]))
LIBS="$XML_LIBS $LIBS"

AX_CHECK_ZLIB

AX_BOOST_BASE([1.34],,AC_MSG_ERROR([suitable Boost not found]))
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"



# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])
AX_CXX_HEADER_TR1_UNORDERED_MAP
AX_CXX_HEADER_TR1_UNORDERED_SET
PKG_CHECK_MODULES([UFC],[ufc-1 = 1.1])
CPPFLAGS="$CPPFLAGS $UFC_CFLAGS"

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


AC_CONFIG_SRCDIR(src/dolfin.h)
AC_CONFIG_HEADERS([src/config/dolfin_config.h])
AC_CONFIG_FILES([Makefile src/Makefile\
			  src/common/Makefile\
			  src/elements/Makefile\
			  src/fem/Makefile\
			  src/function/Makefile\
			  src/graph/Makefile\
			  src/io/Makefile\
			  src/la/Makefile\
			  src/log/Makefile\
			  src/main/Makefile\
			  src/math/Makefile\
			  src/mesh/Makefile\
			  src/mf/Makefile\
			  src/nls/Makefile\
			  src/ode/Makefile\
			  src/parameter/Makefile\
			  src/pde/Makefile\
			  src/quadrature/Makefile])

AC_OUTPUT

# Display help text
echo \
"
---------------------------------------------------------

Configuration of DOLFIN ${PACKAGE_VERSION} finished.

Prefix: '${prefix}'.
Compiler: ${CXX} ${CXXFLAGS} ${CPPFLAGS}

Package features:

  MPI: ${enable_mpi}
  OpenMP: ${enable_openmp}

---------------------------------------------------------
"


