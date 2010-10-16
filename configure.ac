
AC_PREREQ(2.59)
AC_INIT([dolfin],[0.8.1-hpc])
AC_DEFINE_UNQUOTED(DOLFIN_VERSION,["${PACKAGE_VERSION}"],[Version of DOLFIN])
AM_INIT_AUTOMAKE($PACKAGE_NAME, $PACKAGE_VERSION, no-define)
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
PACKAGE_REQUIRES=""
AC_SUBST(PACKAGE_REQUIRES)
AM_SILENT_RULES

# Feature options.
AC_ARG_ENABLE(function-cache, AC_HELP_STRING([--enable-function-cache],
			      [Compile with local data cache in functions]),
			      [enable_function_cache=${enableval}],[enable_function_cache=no])

AC_ARG_ENABLE(optimize-p1, AC_HELP_STRING([--enable-optimize-p1],
			   [Compile with optimization for P1 elements]),
			   [enable_optimize_p1=${enableval}],[enable_optimize_p1=no])

AC_ARG_ENABLE(ublas, AC_HELP_STRING([--disable-ublas], [Compile without support for uBlas]),
		     [enable_ublas=no],[enable_ublas=yes])

AC_ARG_ENABLE(progress-bar, AC_HELP_STRING([--disable-progress-bar],[Compile without progress bar]),
			    [enable_progressbar=no],[enable_progressbar=yes])
AC_ARG_ENABLE(mpi,
              AC_HELP_STRING([--enable-mpi],
	      [Compile with support for MPI]),[enable_mpi=${enableval}], [enable_mpi=no])

AC_ARG_ENABLE(openmp,
              AC_HELP_STRING([--enable-openmp],
	      [Compile with support for OpenMP]),[enable_openmp=${enableval}], [enable_openmp=no])

AC_ARG_ENABLE(janpack, AC_HELP_STRING([--enable-janpack], [Compile with support for JANPACK]),
		     [enable_janpack=${enableval}],[enable_janpack=no])

AC_ARG_ENABLE(petsc, AC_HELP_STRING([--enable-petsc], [Compile with support for PETSc linear algebra]),
		     [enable_petsc=${enableval}],[enable_petsc=no])


if test "x${enable_function_cache}" = xyes; then
   AC_DEFINE(ENABLE_FUNCTION_CACHE, [1], [Compile with local data cache in functions])
fi

if test "x${enable_optimize-p1}" = xyes; then
   AC_DEFINE(ENABLE_P1_OPTIMIZATIONS, [1], [Compile with optimization for P1 elements])
fi

if test "x${enable_ublas}" != xyes; then
   AC_DEFINE(NO_UBLAS, [1], [Compile without support for uBlas])
fi

if test "x${enable_progressbar}" != xyes; then
   AC_DEFINE(NO_PROGRESS_BAR, [1], [Compile without progress bar])
fi



# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
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

if test "x${enable_petsc}" = xyes; then
   PKG_CHECK_MODULES([PETSC],[petsc >= 2.3.3])
   CPPFLAGS="$CPPFLAGS $PETSC_CFLAGS"
   LIBS="$LIBS $PETSC_LIBS"
   AC_DEFINE(HAVE_PETSC,[1],[Have PETSc la-backend])
   AC_SUBST(HAVE_PETSC)
fi

if test "x${enable_janpack}" = xyes; then
   AX_JANPACK
fi

if test "x${enable_mpi}" = xyes; then
   AX_PARMETIS
   LIBS="$PARMETIS_LIBS $LIBS"
fi




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
			  src/quadrature/Makefile\
			  src/dolfin.pc])

AC_OUTPUT

# Display help text
echo \
"
---------------------------------------------------------

Configuration of DOLFIN ${PACKAGE_VERSION} finished.

Prefix: '${prefix}'.
Compiler: ${CXX} ${CXXFLAGS} ${CPPFLAGS}
Libs: ${LD} ${LDFLAGS} ${LIBS}

Package features:

  MPI: ${enable_mpi}
  OpenMP: ${enable_openmp}
  uBlas: ${enable_ublas}
  PETSc: ${enable_petsc}
  JANPACK: ${enable_janpack}
  Function cache: ${enable_function_cache}
  Optimize P1: ${enable_optimize_p1}
  Progress bar: ${enable_progressbar}


---------------------------------------------------------
"


