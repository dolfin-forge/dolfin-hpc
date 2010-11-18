#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_JANPACK],[
	PKG_CHECK_MODULES([janpack], [janpack])
	CPPFLAGS="$CPPFLAGS $janpack_CFLAGS"
	LIBS="$LIBS $janpack_LIBS"
	AC_DEFINE(HAVE_JANPACK,[1],[Define if you have the JANPACK library.])
])