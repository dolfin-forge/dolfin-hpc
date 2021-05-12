
AC_PREREQ(2.59)
AC_INIT([dolfin],[0.9.5-hpc])
AC_DEFINE_UNQUOTED(DOLFIN_VERSION,["${PACKAGE_VERSION}"],[Version of DOLFIN])
AM_INIT_AUTOMAKE([subdir-objects])
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
PACKAGE_REQUIRES=""
AC_SUBST(PACKAGE_REQUIRES)

# Feature options.
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],
             [Compile with debug symbols]),
             [enable_debug=${enableval}],[enable_debug=no])

AC_ARG_ENABLE(function-cache, AC_HELP_STRING([--enable-function-cache],
              [Compile with local data cache in functions]),
              [enable_function_cache=${enableval}],
              [enable_function_cache=no])

AC_ARG_ENABLE(optimize-p1, AC_HELP_STRING([--enable-optimize-p1],
              [Compile with optimization for P1 elements]),
              [enable_optimize_p1=${enableval}],
              [enable_optimize_p1=no])

AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi],
              [Compile with support for MPI]),
              [enable_mpi=${enableval}], [enable_mpi=no])

AC_ARG_ENABLE(mpi-io, AC_HELP_STRING([--enable-mpi-io],
              [Compile with support for MPI I/O]),
              [enable_mpiio=${enableval}],[enable_mpiio=no])

AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp],
              [Compile with support for OpenMP]),
              [enable_openmp=${enableval}], [enable_openmp=no])

AC_ARG_ENABLE(quadrature, AC_HELP_STRING([--enable-quadrature],
              [Enable various quadrature rules]),
              [enable_quad=${enableval}],[enable_quad=no])

AC_ARG_ENABLE(check, AC_HELP_STRING([--enable-check],
              [Compile with support for unit testing, requires libcheck]),
              [enable_check=${enableval}],[enable_check=no])

AC_ARG_ENABLE(python, AC_HELP_STRING([--enable-python],
              [Install Python modules, requires a Python interpreter]),
              [enable_python=${enableval}],[enable_python=no])

AC_ARG_ENABLE(opt-basis, AC_HELP_STRING([--enable-opt-basis],
              [Enable basis function optimization]),
              [enable_opt_basis=${enableval}], [enable_opt_basis=no])

# Optional packages
AC_ARG_WITH(gts, AC_HELP_STRING([--with-gts], [Compile with support for GTS]),
            [with_gts=${withval}],[with_gts=no])

AC_ARG_WITH(janpack, AC_HELP_STRING([--with-janpack],
             [Compile with support for JANPACK]),
             [with_janpack=${withval}],[with_janpack=no])

AC_ARG_WITH(petsc, AC_HELP_STRING([--with-petsc],
           [Compile with support for PETSc linear algebra]),
           [with_petsc=${withval}],[with_petsc=no])

AC_ARG_WITH(trilinos, AC_HELP_STRING([--with-trilinos],
           [Compile with support for Trilinos linear algebra]),
           [with_trilinos=${withval}],[with_trilinos=no])

AC_ARG_WITH(libsim, AC_HELP_STRING([--with-libsim],
            [Compile with support for VisIt/libsim]),
            [with_libsim=${withval}],[with_libsim=no])

AC_ARG_WITH(parmetis, AC_HELP_STRING([--with-parmetis],
         [Compile with support for ParMETIS]),
         [with_parmetis=${withval}],[with_parmetis=no])

AC_ARG_WITH(phmap, AC_HELP_STRING([--with-phmap],
         [Compile with support for parallel hashmap (https://git.io/JfE2p)]),
         [with_phmap=${withval}],[with_phmap=no])

AC_ARG_WITH(zoltan, AC_HELP_STRING([--with-zoltan],
         [Compile with support for Zoltan]),
         [with_zoltan=${withval}],[with_zoltan=no])

AC_ARG_WITH(convert, AC_HELP_STRING([--with-convert],
         [Build external mesh converter tool]),
         [with_convert=${withval}],[with_convert=no])


if test "x${enable_function_cache}" = xyes; then
   AC_DEFINE(ENABLE_FUNCTION_CACHE, [1],
   [Compile with local data cache in functions])
fi

if test "x${enable_optimize_p1}" = xyes; then
   AC_DEFINE(ENABLE_P1_OPTIMIZATIONS, [1],
   [Compile with optimization for P1 elements])
fi

if test "x${enable_opt_basis}" = xyes; then
   AC_DEFINE(ENABLE_EVALUATE_BASIS_FROM_COORDINATES, [1],
   [Enable basis function optimization])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# Ensure strict C++14 compliance
AX_CXX_COMPILE_STDCXX_14([noext], [mandatory])

# add default warning flags
CXXFLAGS="-Wall -Wextra -pedantic $CXXFLAGS"

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
  CXXFLAGS="$OPENMP_CXXFLAGS $CXXFLAGS"
else
   if test "x${ax_cv_cxx_compiler_vendor}" = xcray; then
      CXXFLAGS="-hnoomp $CXXFLAGS"
   fi
fi

# Store build information (date, host, CXX)
AC_DEFINE_UNQUOTED(DOLFIN_BUILD_INFO,
  ["(Build: $(date +%Y-%m-%d) on $host using $ax_cv_cxx_compiler_vendor $ax_cv_cxx_compiler_version)"], [])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# UFL is now default
PKG_CHECK_MODULES([UFC],[ufc-1 >= 2.4.0])
CPPFLAGS="$CPPFLAGS $UFC_CFLAGS"

# Checks for libraries.
AX_ZLIB

if test "${ZLIB_LIBS}"; then
   LIBS="$ZLIB_LIBS $LIBS"
   enable_zlib="yes"
else
   enable_zlib="no"
fi

AX_CRAY

if test "x${enable_mpi}" != xno; then

found_domain_decomp="no"

   # This assumes a modern CRAY, for example XT, XE or XC series systems
   if test "x${is_cray}" = xyes; then
      if test "x${with_parmetis}" != xno; then
         AX_CRAY_PARMETIS
   if test "x${have_cray_parmetis}" = xyes; then
      found_domain_decomp="yes"
   else
      with_parmetis="no"
   fi
      fi

      if test "x${with_zoltan}" != xno; then
         if test "x${enable_optimize_p1}" != xno; then
      AC_MSG_ERROR([Zoltan requires facet connectivity])
   fi

         AX_CRAY_ZOLTAN
   if test "x${have_cray_zoltan}" = xyes; then
      found_domain_decomp="yes"
   else
      with_zoltan="no"
   fi
      fi
   else
      if test "x${with_parmetis}" != xno; then
         AX_PARMETIS
         if test "${PARMETIS_LIBS}"; then
            LIBS="$PARMETIS_LIBS $LIBS"
      found_domain_decomp="yes"
      with_parmetis="yes"
   else
      with_parmetis="no"
         fi
      fi

      if test "x${with_zoltan}" != xno; then
         if test "x${enable_optimize_p1}" != xno; then
      AC_MSG_ERROR([Zoltan requires facet connectivity])
   fi

         AX_ZOLTAN
         if test "${ZOLTAN_LIBS}"; then
            LIBS="$ZOLTAN_LIBS $LIBS"
      found_domain_decomp="yes"
      with_zoltan="yes"
   else
      with_zoltan="no"
         fi
      fi
   fi

   if test "x${found_domain_decomp}" != xyes; then
      AC_MSG_ERROR([Required partitioning library not found])
   fi
fi

if test "x${with_gts}" = xyes; then
   PKG_CHECK_MODULES([gts], [gts])
   CPPFLAGS="$CPPFLAGS $gts_CFLAGS"
   LIBS="$LIBS $gts_LIBS"
   AC_DEFINE(HAVE_GTS,[1],[Have GTS library])
fi

if test "x${with_phmap}" != xno; then
  AX_PHMAP
fi

found_la_backend="no"

if test "x${with_petsc}" != xno; then
   # This assumes a modern CRAY, for example XT, XE or XC series systems
   if test "x${is_cray}" = xyes; then
      AX_CRAY_PETSC
      if test "x${have_cray_petsc}" = xyes; then
         found_la_backend="yes"
         with_petsc="yes"
      else
   with_petsc="no"
      fi
   else
      AX_PETSC
      if test "${PETSC_LDFLAGS}"; then
         found_la_backend="yes"
         with_petsc="yes"
      else
         with_petsc="no"
      fi
   fi
fi

if test "x${with_trilinos}" != xno; then
  AX_TRILINOS
fi

if test "x${with_janpack}" != xno; then
   AX_JANPACK
   if test "${janpack_LIBS}"; then
     found_la_backend="yes"
   elif test "${janpack_mpi_LIBS}"; then
     found_la_backend="yes"
   else
     with_janpack="no"
   fi
fi

if test "x${found_la_backend}" != xyes; then
   AC_MSG_ERROR([Linear Algebra backend not found])
fi

if test "x${enable_quad}" = xyes; then
   # This assumes a modern CRAY, for example XT or XE series systems
   if test "x${ax_cv_cxx_compiler_vendor}" = xcray; then
      AX_CRAY_LIBSCI
      if test $have_cray_libsci = yes; then
         have_lapack=yes
      else
   have_lapack=no
      fi
   elif test "x${ax_cv_cxx_compiler_vendor}" = xsgi; then
      AC_CHECK_LIB(scs,dgesv_,[have_lapack=yes, LIBS="$LIBS -lscs"],
                   have_lapack=no,)
   elif test "x${ax_cv_cxx_compiler_vendor}" = xsun; then
      AC_CHECK_LIB(sunperf,dgesv_,[have_lapack=yes, LIBS="$LIBS -lsunperf"],
                   [have_lapack=no],)
   else
      AC_CHECK_LIB(lapack,dgesv_,
                   [have_lapack=yes, LIBS="$LIBS -lblas -llapack"],
                   [have_lapack=no], [-lblas])
   fi

   if test $have_lapack = no; then
      AC_CHECK_LIB(lapack,dgesv_,[have_lapack=yes, LIBS="$LIBS"],
                   have_lapack=no,)
   fi


   if test $have_lapack = yes; then
      AC_DEFINE(HAVE_F77_LAPACK, [1], [Have FORTRAN LAPACK])
   elif test $have_lapack = no; then
      AC_MSG_ERROR([LAPACK not found, needed for various quadrature rules])
   fi
fi

if test "x${with_libsim}" != xno; then
   AX_LIBSIM
fi

# Test for unit testing framework
if test "x${enable_check}" = xyes; then
   PKG_CHECK_EXISTS([check],[have_check=yes],[have_check=no])
   if test "x${have_check}" = xyes; then
      PKG_CHECK_MODULES([CHECK], [check >= 0.12.0])
      AC_DEFINE(HAVE_CHECK,[1],[Have CHECK library])
   else
      AC_MSG_ERROR([CHECK not found, needed for unit testing])
   fi
fi
# Add variable for automake to compile tests
AM_CONDITIONAL([ENABLE_CHECK], [test "x${enable_check}" = xyes])

if test "x${enable_debug}" != xno; then
   CXXFLAGS="-g -DDEBUG=1 $CXXFLAGS"
   case "$host_os" in
   *solaris*)
       LIBS="$LIBS -lumem"
       AC_MSG_RESULT([[Using DEBUG mode on Illumos/Solaris with libumem support]])
       ;;
   *)
   AC_MSG_RESULT([[Using DEBUG mode]])
       ;;
   esac
else
   AC_MSG_RESULT([[Using RELEASE mode]])
   CXXFLAGS="-DNDEBUG $CXXFLAGS"
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

# Checks for Python
if test "x${enable_python}" = xyes; then
    AM_PATH_PYTHON([3.0])
    if ! ${HAVE_PYTHON}; then
        AC_MSG_ERROR([Python interpreter not found, module cannot be installed.])
    fi
    case ${PYTHON_VERSION} in
    3*)
    ;;
    *)
    AC_MSG_ERROR([No Python3 interpreter found, module cannot be installed.])
    ;;
    esac
fi

# Add variable for installation of Python modules
AM_CONDITIONAL([ENABLE_PYTHON_MODULES], [test "x${enable_python}" = xyes])

# Checks for Doxygen
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -n "$DOXYGEN"; then
  AC_CHECK_PROGS([DOXYGEN_DEPS], [dot])
fi
AM_CONDITIONAL([ENABLE_DOXYGEN], [test -n "$DOXYGEN_DEPS"])

AM_CONDITIONAL([ENABLE_MESH_CONVERTER], false)
# Configure external mesh converting tool
if test "x${with_convert}" != xno; then
   AC_CONFIG_SUBDIRS(./misc/utils/convert/)
   AM_CONDITIONAL([ENABLE_MESH_CONVERTER], true)
fi

AC_CONFIG_HEADERS([include/dolfin/config/dolfin_config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile\
                 src/adaptivity/Makefile\
                 src/common/Makefile\
                 src/evolution/Makefile\
                 src/fem/Makefile\
                 src/function/Makefile\
                 src/io/Makefile\
                 src/insitu/Makefile\
                 src/la/Makefile\
                 src/la/janpack/Makefile\
                 src/la/petsc/Makefile\
                 src/la/trilinos/Makefile\
                 src/log/Makefile\
                 src/main/Makefile\
                 src/math/Makefile\
                 src/mesh/Makefile\
                 src/mesh/celltypes/Makefile\
                 src/mesh/entities/Makefile\
                 src/mesh/mappings/Makefile\
                 src/mesh/unitmeshes/Makefile\
                 src/mesh/utilities/Makefile\
                 src/parameter/Makefile\
                 src/quadrature/Makefile\
                 src/dolfin.pc\
                 include/Makefile\
                 include/dolfin/Makefile])

# Tests
AC_CONFIG_FILES([tests/Makefile\
                 tests/common/Makefile\
                 tests/fem/Makefile\
                 tests/function/Makefile\
                 tests/io/Makefile\
                 tests/la/Makefile\
                 tests/log/Makefile\
                 tests/main/Makefile\
                 tests/mesh/Makefile\
                 tests/parameter/Makefile])

# Python
AC_CONFIG_FILES([site-packages/Makefile])

# Doxygen
AC_CONFIG_FILES([doc/progr_reference/Doxyfile\
                 doc/progr_reference/Makefile])

AC_OUTPUT

# Display help text
echo \
"
---------------------------------------------------------

Configuration of DOLFIN ${PACKAGE_VERSION} finished.

Prefix: '${prefix}'.
Compiler: ${CXX} ${CPPFLAGS} ${CXXFLAGS}
Compiler: ${CXX}
Compiler: ${CXXFLAGS}
Compiler: ${CPPFLAGS}
Libs: ${LD} ${LDFLAGS} ${LIBS}

Package features:

  MPI:                         ${enable_mpi}
  MPI I/O:                     ${enable_mpiio}
  OpenMP:                      ${enable_openmp}
  PETSc:                       ${with_petsc}
  Trilinos:                    ${with_trilinos}
  JANPACK:                     ${with_janpack}
  GTS:                         ${with_gts}
  ParMETIS:                    ${with_parmetis}
  Zoltan:                      ${with_zoltan}
  libsim:                      ${with_libsim}
  Function cache:              ${enable_function_cache}
  Optimize P1:                 ${enable_optimize_p1}
  Compressed VTK:              ${enable_zlib}
  Quadrature:                  ${enable_quad}
  Parallel Hash Map:           ${have_phmap_h}
  Basis function optimization: ${enable_opt_basis}
  CHECK support:               ${enable_check}
  DEBUG mode:                  ${enable_debug}

---------------------------------------------------------
"
if test "x${enable_python}" = xyes; then
echo \
"  Python modules will be installed in:

        ${pythondir}

  with \${prefix} = ${prefix} .

  Do not forget to set the PYTHONPATH variable.
---------------------------------------------------------
"
fi


