#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_ZLIB],[
	AC_ARG_WITH([zlib],
	AS_HELP_STRING([--with-zlib=DIR],
	[Directory for zlib]),
	[	   
	if test -d "$withval"; then
		ac_zlib_path="$withval";
		ZLIB_LDFLAGS="-L$ac_zlib_path/lib"  
		ZLIB_CPPFLAGS="-I$ac_zlib_path/include"
	fi
	],)

	AC_ARG_WITH([zlib-libdir],
	AS_HELP_STRING([--with-zlib-libdir=LIBDIR],
	[Directory for zlib library]),
	[
	if test -d "$withval"; then
	   ac_zlib_libdir="$withval"
	   ZLIB_LDFLAGS="-L$ac_zlib_libdir"  	   
	fi
	],)

	if test -d "$ac_zlib_path"; then
	   CPPFLAGS_SAVED="$CPPFLAGS"
	   LDFLAGS_SAVED="$LDFLAGS"
	   CPPFLAGS="$ZLIB_CPPFLAGS $CPPFLAGS"
	   LDFLAGS="$ZLIB_LDFLAGS $LDFLAGS"
	   export CPPFLAGS
	   export LDFLAGS
	fi

	if test -d "$ac_zlib_libdir"; then	   
	    LDFLAGS="$ZLIB_LDFLAGS $LDFLAGS"
	fi

			
	AC_CHECK_HEADER([zlib.h],[have_zlib_h=yes],[have_zlib_h=no])
	AC_CHECK_LIB(z, inflateEnd,[have_zlib=yes;ZLIB_LIBS="-lz"],[have_zlib=no],)
	AC_SUBST(ZLIB_LIBS)
	if test x"${have_zlib}" = xyes; then
	   AC_DEFINE(HAVE_LIBZ,1,[Define if you have the Zlib library.])
	else
		if test -d "$ac_zlib_path"; then	
		   CPPFLAGS="$CPPFLAGS_SAVED"
		   LDFLAGS="$LDFLAGS_SAVED"
		fi
	fi
])


