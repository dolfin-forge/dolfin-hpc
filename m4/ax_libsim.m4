#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <leifniclas.jansson@riken.jp> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet
# some day, and you think this stuff is worth it, you can buy me a
# beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_LIBSIM],[
	AC_ARG_WITH(libsim,
	AC_HELP_STRING([--with-libsim],
	      [Compile with support for VisIt/libsim in-situ visualization]),
	[	      	     
      	if test -d "$withval"; then
		ac_libsim_dir="$withval"		
		LIBSIM_CPPFLAGS="-I$ac_libsim_dir/V2/include"
		LIBSIM_LDFLAGS="-L$ac_libsim_dir/V2/lib"
		CPPFLAGS_SAVED="$CPPFLAGS"
		LDFLAGS_SAVED="$LDFLAGS"
		CPPFLAGS="$LIBSIM_CPPFLAGS $CPPFLAGS"
		LDFLAGS="$LIBSIM_LDFLAGS $LDFLAGS"
		export CPPFLAGS
		export LDFLAGS
		with_libsim=yes
	fi
	],[with_libsim=no])

	if test "x$with_libsim" != xno; then
	   AC_MSG_CHECKING(for VisIt/libsim)
	   AC_COMPILE_IFELSE([
	   #include <VisItControlInterface_V2.h>	  
	   void main() {
	   VisItSetupEnvironment2(0);
	   }
	   ],[have_libsim=yes], [have_libsim=no])

          if test "x$have_libsim" = "xyes"; then
   	      AC_MSG_RESULT([yes])
	      LIBS="$LIBS -lsimV2 -ldl"
	      AC_DEFINE([HAVE_LIBSIM], [1], [VisIt/Libsim support])
	  else
	      AC_MSG_RESULT([no])
              CPPFLAGS="$CPPFLAGS_SAVED"
              LDFLAGS="$LDFLAGS_SAVED"

	  fi
        fi	  
       ])