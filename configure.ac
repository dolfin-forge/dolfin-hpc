
AC_PREREQ(2.59)
AC_INIT([dolfin],[0.8.1-hpc])
AC_DEFINE_UNQUOTED(DOLFIN_VERSION,["${PACKAGE_VERSION}"],[Version of DOLFIN])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
PACKAGE_REQUIRES=""
AC_SUBST(PACKAGE_REQUIRES)



# Feature options.
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[Compile with debug symbols]),
		     [enable_debug=${enableval}],[enable_debug=no])

AC_ARG_ENABLE(function-cache, AC_HELP_STRING([--enable-function-cache],
			      [Compile with local data cache in functions]),
			      [enable_function_cache=${enableval}],[enable_function_cache=no])

AC_ARG_ENABLE(optimize-p1, AC_HELP_STRING([--enable-optimize-p1],
			   [Compile with optimization for P1 elements]),
			   [enable_optimize_p1=${enableval}],[enable_optimize_p1=no])

AC_ARG_ENABLE(progress-bar, AC_HELP_STRING([--disable-progress-bar],[Compile without progress bar]),
			    [enable_progressbar=no],[enable_progressbar=yes])
AC_ARG_ENABLE(mpi,
              AC_HELP_STRING([--enable-mpi],
	      [Compile with support for MPI]),[enable_mpi=${enableval}], [enable_mpi=no])

AC_ARG_ENABLE(mpi-io, AC_HELP_STRING([--enable-mpi-io],
		     [Compile with support for MPI I/O]),
		     [enable_mpiio=yes],[enable_mpiio=no])

AC_ARG_ENABLE(openmp,
              AC_HELP_STRING([--enable-openmp],
	      [Compile with support for OpenMP]),[enable_openmp=${enableval}], [enable_openmp=no])

AC_ARG_ENABLE(boost-tr1, AC_HELP_STRING([--enable-boost-tr1],
			 [Compile with support for Boost TR1 headers]),
			 [enable_boosttr1=yes],[enable_boostr1=no])





# Optional packages
AC_ARG_WITH(gts, AC_HELP_STRING([--with-gts], [Compile with support for GTS]),
		 [with_gts=${withval}],[with_gts=no])

AC_ARG_WITH(janpack, AC_HELP_STRING([--with-janpack], [Compile with support for JANPACK]),
		     [with_janpack=${withval}],[with_janpack=no])

AC_ARG_WITH(petsc, AC_HELP_STRING([--with-petsc], [Compile with support for PETSc linear algebra]),
		     [with_petsc=${withval}],[with_petsc=no])


if test "x${enable_function_cache}" = xyes; then
   AC_DEFINE(ENABLE_FUNCTION_CACHE, [1], [Compile with local data cache in functions])
fi

if test "x${enable_optimize_p1}" = xyes; then
   AC_DEFINE(ENABLE_P1_OPTIMIZATIONS, [1], [Compile with optimization for P1 elements])
fi

if test "x${enable_progressbar}" != xyes; then
   AC_DEFINE(NO_PROGRESS_BAR, [1], [Compile without progress bar])
fi

if test "x${enable_boosttr1}" = xyes; then
   AC_DEFINE(ENABLE_BOOST_TR1,[1],[Compile with support for Boost TR1 headers])
   need_boost=yes	
fi


# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG
AX_COMPILER_VENDOR
if test "x${enable_mpi}" != xno; then
   AX_MPI	
   CXX="$MPICXX"
   LIBS="$MPILIBS $LIBS"   
   if test "x${ax_cv_cxx_compiler_vendor}" = xsgi; then
      LIBS="$LIBS -lmpi++"
   fi   

   if test "x${enable_mpiio}" != xno; then
      AC_DEFINE(ENABLE_MPIIO, 1, [Compile with support for MPI I/O])
   fi
fi



if test "x${enable_openmp}" != xno; then
   AX_OPENMP
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
else
   if test "x${ax_cv_cxx_compiler_vendor}" = xcray; then
      CXXFLAGS="-hnoomp $CXXFLAGS" 
   fi
fi



# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])
AX_CXX_HEADER_TR1_UNORDERED_MAP
AX_CXX_HEADER_TR1_UNORDERED_SET
PKG_CHECK_MODULES([UFC],[ufc-1 = 1.1])
CPPFLAGS="$CPPFLAGS $UFC_CFLAGS"


if test "x${need_boost}" = xyes; then
   AX_BOOST_BASE([1.34],,AC_MSG_ERROR([suitable Boost not found]))
   CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
fi

if test "x${enable_boosttr1}" = xyes; then
   if test "x${ax_cv_cxx_compiler_vendor}" = xsun; then
      CXXFLAGS="-D_RWSTD_ALLOCATOR $CXXFLAGS"
   fi
fi


# Checks for libraries.
AM_PATH_XML2(,,AC_MSG_ERROR([libxml2 not found]))
LIBS="$XML_LIBS $LIBS"

AX_ZLIB

if test "${ZLIB_LIBS}"; then
   LIBS="$ZLIB_LIBS $LIBS"
   enable_zlib="yes"
else
   enable_zlib="no"
fi


if test "x${enable_mpi}" != xno; then
   AX_PARMETIS
   if test "${PARMETIS_LIBS}"; then
      LIBS="$PARMETIS_LIBS $LIBS"
   else
      AC_MSG_ERROR([Required library ParMetis not found])
   fi
fi

if test "x${with_gts}" = xyes; then
   PKG_CHECK_MODULES([gts], [gts])
   CPPFLAGS="$CPPFLAGS $gts_CFLAGS"
   LIBS="$LIBS $gts_LIBS"
   AC_DEFINE(HAVE_GTS,[1],[Have GTS library])
fi

if test "x${with_petsc}" != xno; then
   # This assumes a modern CRAY, for example XT or XE series systems
   if test "x${ax_cv_cxx_compiler_vendor}" = xcray; then
      AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
   else
      AX_PETSC
   fi
fi

if test "x${with_janpack}" != xno; then
   AX_JANPACK
fi

if test "x${enable_debug}" != xno; then
   CXXFLAGS="-g -DDEBUG $CXXFLAGS"
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

if test "x${ac_cv_c_bigendian}" != xno; then
   AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Big endian byte order.])
fi


# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor memset pow sqrt])


AC_CONFIG_HEADERS([include/dolfin/config/dolfin_config.h])
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
			  src/nls/Makefile\
			  src/parameter/Makefile\
			  src/pde/Makefile\
			  src/dolfin.pc
			  include/Makefile
			  include/dolfin/Makefile])

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
  MPI I/O: ${enable_mpiio}
  OpenMP: ${enable_openmp}
  PETSc: ${with_petsc}
  JANPACK: ${with_janpack}
  GTS: ${with_gts}
  Function cache: ${enable_function_cache}
  Optimize P1: ${enable_optimize_p1}
  Progress bar: ${enable_progressbar}
  Compressed VTK: ${enable_zlib}


---------------------------------------------------------
"


