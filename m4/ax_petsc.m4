#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_PETSC],[
	AC_MSG_CHECKING(for PETSc)
	AC_LANG(C++)
	AC_ARG_WITH([petsc-dir],
	AS_HELP_STRING([--with-petsc-dir=DIR],
	[Directory for petsc]),
	[	   
	if test -d "$withval"; then
		ac_petsc_dir="$withval"
		PETSC_DIR="$withval"
		export PETSC_DIR
		AC_SUBST(PETSC_DIR)
	fi
	],)


	CPPFLAGS_SAVED="$CPPFLAGS"
	LDFLAGS_SAVED="$LDFLAGS"

	if test -d "$ac_petsc_dir"; then
	   if test -d "$ac_petsc_dir/bmake"; then
	      ac_petsc_arch=`grep PETSC_ARCH  $ac_petsc_dir/bmake/petscconf | sed 's/PETSC_ARCH=/''/'`
	      ac_petsc_libdir="$ac_petsc_dir/lib/$ac_petsc_arch"	      	      
	      PETSC_CPPFLAGS="-I$ac_petsc_dir/bmake/$ac_petsc_arch -I$ac_petsc_dir/include"
	      PETSC_LDFLAGS="-L$ac_petsc_libdir -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc"
	      
	      have_petsc="yes"
	   elif test -d "$ac_petsc_dir/conf"; then
	     ac_petsc_libdir="$ac_petsc_dir/lib"
	     PETSC_PKG_LIBS=`grep PACKAGES_LIBS $ac_petsc_dir/conf/petscvariables | sed 's/PACKAGES_LIBS =/''/'`
	     PETSC_CPPFLAGS="-I$ac_petsc_dir/include"
	     PETSC_LDFLAGS="-L$ac_petsc_libdir $PETSC_PKG_LIBS -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc"

	     have_petsc="yes"
	   fi
		
	fi

	if test x"${have_petsc}" = xyes; then
	   AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
	   CPPFLAGS="$CPPFLAGS $PETSC_CPPFLAGS"
	   LDFLAGS="$LD_FLAGS $PETSC_LDFLAGS"
	   AC_MSG_RESULT(yes)
	else
	   CPPFLAGS="$CPPFLAGS_SAVED"
	   LDFLAGS="$LDFLAGS_SAVED"
           AC_MSG_RESULT(no)
	fi
])


