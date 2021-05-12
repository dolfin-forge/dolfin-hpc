#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <julian@ingridcloud.com> wrote this file. As long as you retain this notice
# you can do whatever you want with this stuff. If we meet some day, and you
# think this stuff is worth it, you can buy me a beer in return Julian Hornich
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_TRILINOS],[
	AC_MSG_CHECKING(for Trilinos)
	AC_ARG_WITH([trilinos],
	AS_HELP_STRING([--with-trilinos=DIR],
	[Directory for trilinos]),
	[
	if test -d "$withval"; then
		ac_trilinos_path="$withval";
	else
	  ac_trilinos_path="/usr"
	fi
	],)

	if test -f "${ac_trilinos_path}/include/Trilinos_version.h"; then
	  TRILINOS_MAKEFILE=${ac_trilinos_path}/include/Makefile.export.Trilinos

    # add trilinos CXXFLAGS and include dirs
	  CXXFLAGS="$(grep Trilinos_INCLUDE_DIRS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $CXXFLAGS"
	  CXXFLAGS="$(grep Trilinos_TPL_INCLUDE_DIRS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $CXXFLAGS"
	  CXXFLAGS="$(grep Trilinos_CXX_COMPILER_FLAGS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $CXXFLAGS"

    # add trilinos link flags and dirs
	  LDFLAGS="$(grep Trilinos_LIBRARY_DIRS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $LDFLAGS"
	  LDFLAGS="$(grep Trilinos_TPL_LIBRARY_DIRS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $LDFLAGS"
	  LDFLAGS="$(grep Trilinos_EXTRA_LD_FLAGS ${TRILINOS_MAKEFILE} | sed 's/.*=//') $LDFLAGS"

    # add trilinos libraries
	  LIBS="$(grep Trilinos_TPL_LIBRARIES ${TRILINOS_MAKEFILE} | sed 's/.*=//') $LIBS"
	  LIBS="$(grep Trilinos_LIBRARIES ${TRILINOS_MAKEFILE} | sed 's/.*=//') $LIBS"

    AC_DEFINE(HAVE_TRILINOS, 1, [Define if you have the trilinos library.])
    found_la_backend="yes"

    AC_MSG_RESULT(yes (version $(grep Trilinos_VERSION ${TRILINOS_MAKEFILE} | sed 's/.*= *//')))
    with_trilinos=yes
  else
    AC_MSG_RESULT(no)
    with_trilinos=no
	fi
])
