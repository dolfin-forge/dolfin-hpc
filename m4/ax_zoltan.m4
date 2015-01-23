#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_ZOLTAN],[
	AC_ARG_WITH([zoltan],
	AS_HELP_STRING([--with-zoltan=DIR],
	[Directory for zoltan]),
	[	   
	if test -d "$withval"; then
		ac_zoltan_path="$withval";
		ZOLTAN_LDFLAGS="-L$ac_zoltan_path/lib"  
		ZOLTAN_CPPFLAGS="-I$ac_zoltan_path/include"
	fi
	],)

	AC_ARG_WITH([zoltan-libdir],
	AS_HELP_STRING([--with-zoltan-libdir=LIBDIR],
	[Directory for zoltan library]),
	[
	if test -d "$withval"; then
	   ac_zoltan_libdir="$withval"
	fi
	],)

	if test -d "$ac_zoltan_libdir"; then	   
	    ZOLTAN_LDFLAGS="-L$ac_zoltan_libdir"  	   
        fi

	if test -d "$ac_zoltan_path"; then
	   CPPFLAGS_SAVED="$CPPFLAGS"
	   LDFLAGS_SAVED="$LDFLAGS"
	   CPPFLAGS="$ZOLTAN_CPPFLAGS $CPPFLAGS"
	   LDFLAGS="$ZOLTAN_LDFLAGS $LDFLAGS"
	   export CPPFLAGS
	   export LDFLAGS
	fi
			
	AC_CHECK_HEADER([zoltan.h],[have_zoltan_h=yes],[have_zoltan_h=no])
	AC_CHECK_LIB(zoltan, Zoltan_LB_Partition, 
			     [have_zoltan=yes;ZOLTAN_LIBS="-lzoltan"],
			     [have_zoltan=no])
	AC_SUBST(ZOLTAN_LIBS)
	if test x"${have_zoltan}" = xyes; then
	   AC_DEFINE(HAVE_ZOLTAN,1,[Define if you have the Zoltan library.])
	else
		if test -d "$ac_zoltan_path"; then	
		   CPPFLAGS="$CPPFLAGS_SAVED"
		   LDFLAGS="$LDFLAGS_SAVED"
		fi
	fi
])


