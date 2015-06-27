#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_CRAY],[
	AC_MSG_CHECKING([for a Cray XT, XE, XC system])
	AC_EGREP_CPP(yes,
	[#if defined(__CRAYXT) || defined(__CRAYXE) || defined(__CRAYXC)
	  yes
	 #endif
	],
	[AC_MSG_RESULT([yes])
	is_cray="yes"],
	[is_cray="no"
	AC_MSG_RESULT([no])])
	AC_SUBST(is_cray)])

AC_DEFUN([AX_CRAY_PETSC],[
	AC_MSG_CHECKING([Cray PETSc])
	if test "${CRAY_PETSC_VERSION}"; then
	   have_cray_petsc="yes"
	else
	   have_cray_petsc="no"
	fi
	AC_SUBST(have_cray_petsc)
	if test "x${have_cray_petsc}" = xyes; then
	   AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
	   AC_MSG_RESULT([yes])
	else
	   AC_MSG_RESULT([no])
	fi
])

AC_DEFUN([AX_CRAY_PARMETIS],[
	AC_MSG_CHECKING([Cray ParMETIS])
	if test "${CRAY_TRILINOS_VERSION}"; then
	   have_cray_parmetis="yes"
	elif test "${CRAY_PETSC_VERSION}"; then
	   have_cray_parmetis="yes"
	else
	   have_cray_parmetis="no"
	fi
	AC_SUBST(have_cray_parmetis)
	if test "x${have_cray_parmetis}" = xyes; then
           AC_DEFINE(HAVE_PARMETIS,1,
		     [Define if you have the ParMETIS library])
	   AC_MSG_RESULT([yes])
	else
	   AC_MSG_RESULT([no])
	fi
])

AC_DEFUN([AX_CRAY_ZOLTAN],[
	AC_MSG_CHECKING([Cray Zoltan])
	if test "${CRAY_TRILINOS_VERSION}"; then
	   have_cray_zoltan="yes"
	else
	   have_cray_zoltan="no"
	fi
	AC_SUBST(have_cray_zoltan)
	if test "x${have_cray_zoltan}" = xyes; then
	   AC_DEFINE(HAVE_ZOLTAN,1,
		     [Define if you have the Zoltan library.])
	   AC_MSG_RESULT([yes])
	else
	   AC_MSG_RESULT([no])
	fi
])

