#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_PARMETIS],[
	AC_CHECK_HEADER([parmetis.h],[have_parmetis_h=yes],[have_parmetis_h=no])
	AC_CHECK_LIB(parmetis, ParMETIS_V3_PartMeshKway,[have_parmetis=yes;PARMETIS_LIBS="-lparmetis -lmetis"],[have_parmetis=no],[-lmetis])
	AC_SUBST(PARMETIS_LIBS)
	if test x"${have_parmetis}" = xyes; then
	   AC_DEFINE(HAVE_PARMETIS,1,[Define if you have the ParMETIS library.])
	fi
])


