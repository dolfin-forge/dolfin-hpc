#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_JANPACK],[
	AC_CHECK_LIB(janpack,init_mat_crs,[have_janpack=yes],[have_janpack=no])
	if test x"${have_janpack}" = xyes; then
	   AC_DEFINE(HAVE_JANPACK,1,[Define if you have the JANPACK library.])
	fi
])