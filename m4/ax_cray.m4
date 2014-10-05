#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
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

