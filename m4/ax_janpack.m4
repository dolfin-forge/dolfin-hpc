#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_JANPACK],[
	PKG_CHECK_MODULES([janpack], [janpack >= 0.2.3], 
				     have_janpack=yes, have_janpack=no)
	if test "x${have_janpack}" = xyes; then 	
	   CPPFLAGS="$CPPFLAGS $janpack_CFLAGS"
	   LIBS="$LIBS $janpack_LIBS"
	   AC_DEFINE(HAVE_JANPACK,[1],
		     [Define if you have the JANPACK library.])
	fi

	if test "x${have_janpack}" = xno; then 	
	PKG_CHECK_MODULES([janpack_mpi], [janpack >= mpi-0.2.3], 
					 have_janpack=yes, have_janpack=no)
		if test "x${have_janpack}" = xyes; then 	
		   CPPFLAGS="$CPPFLAGS $janpack_mpi_CFLAGS"
		   LIBS="$LIBS $janpack_mpi_LIBS"
		   AC_DEFINE(HAVE_JANPACK,[1],
			     [Define if you have the JANPACK library.])
		   AC_DEFINE(HAVE_JANPACK_MPI,[1],
			     [Define if you have the MPI version of JANPACK.])
		fi
	fi
])
