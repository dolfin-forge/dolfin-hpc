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
	AC_ARG_WITH([petsc],
	AS_HELP_STRING([--with-petsc=DIR],
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
	have_petsc=no
	have_pkgconfig_petsc="no"
	if test -d "$ac_petsc_dir"; then
	   if test -d "$ac_petsc_dir/bmake"; then
	      ac_petsc_arch=`grep PETSC_ARCH  $ac_petsc_dir/bmake/petscconf | sed 's/PETSC_ARCH=/''/'`
	      ac_petsc_libdir="$ac_petsc_dir/lib/$ac_petsc_arch"	      	      
	      CC_LINKER_SLFLAGS=`grep 'CC_LINKER_SLFLAG ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/CC_LINKER_SLFLAG =/''/'`
	      X11_LIB=`grep 'X11_LIB ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/X11_LIB =/''/'`
	      HYPRE_LIB=`grep 'HYPRE_LIB ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/HYPRE_LIB =/''/'`
	      UMFPACK_LIB=`grep 'UMFPACK_LIB ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/UMFPACK_LIB =/''/'`
	      MPI_LIB=`grep 'MPI_LIB ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/MPI_LIB =/''/'`
	      BLASLAPACK_LIB=`grep 'BLASLAPACK_LIB ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/BLASLAPACK_LIB =/''/'`
	      PETSC_MPI_INC=`grep 'MPI_INCLUDE ='  $ac_petsc_dir/bmake/$ac_petsc_arch/petscconf | sed 's/MPI_INCLUDE =/''/'`
	      PETSC_CPPFLAGS="-I$ac_petsc_dir/bmake/$ac_petsc_arch -I$ac_petsc_dir/include/ $PETSC_MPI_INC"
	      PETSC_LDFLAGS="$CC_LINKER_SLFLAGS$ac_petsc_libdir -L$ac_petsc_libdir -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc  $X11_LIB $HYPRE_LIB $UMFPACK_LIB $MPI_LIB $BLASLAPACK_LIB -lm"
	      
	      have_petsc="yes"
	   elif test -d "$ac_petsc_dir/conf"; then
	   
          cat <<EOF >config_petsc
include $PETSC_DIR/conf/petscvariables

petsclibs:
	echo -L$PETSC_DIR/lib/  \$(PETSC_LIB)
petscinc:
	echo \$(PETSC_CC_INCLUDES)
EOF
	     PETSC_LDFLAGS=`make -s -f config_petsc petsclibs`
	     PETSC_CPPFLAGS=`make -s -f config_petsc petscinc`
	     rm -fr config_petsc
	     have_petsc="yes"
	  elif test -d "$ac_petsc_dir/lib/petsc/conf"; then
	   
          cat <<EOF >config_petsc
include $PETSC_DIR/lib/petsc/conf/petscvariables

petsclibs:
	echo -L$PETSC_DIR/lib/  \$(PETSC_LIB)
petscinc:
	echo \$(PETSC_CC_INCLUDES)
EOF
	     PETSC_LDFLAGS=`make -s -f config_petsc petsclibs`
	     PETSC_CPPFLAGS=`make -s -f config_petsc petscinc`
	     rm -fr config_petsc
	     have_petsc="yes"

	  fi		
	fi
	if test x"${have_petsc}" = xno; then
	   if test -d "$ac_petsc_dir/conf"; then
	   
          cat <<EOF >config_petsc
include $PETSC_DIR/conf/petscvariables

petsclibs:
	echo -L$PETSC_DIR/lib64/  \$(PETSC_LIB)
petscinc:
	echo \$(PETSC_CC_INCLUDES)
EOF
	     PETSC_LDFLAGS=`make -s -f config_petsc petsclibs`
	     PETSC_CPPFLAGS=`make -s -f config_petsc petscinc`
	     rm -fr config_petsc
	     have_petsc="yes"
	  elif test -d "$ac_petsc_dir/lib64/petsc/conf"; then
	   
          cat <<EOF >config_petsc
include $PETSC_DIR/lib/petsc/conf/petscvariables

petsclibs:
	echo -L$PETSC_DIR/lib64/  \$(PETSC_LIB)
petscinc:
	echo \$(PETSC_CC_INCLUDES)
EOF
	     PETSC_LDFLAGS=`make -s -f config_petsc petsclibs`
	     PETSC_CPPFLAGS=`make -s -f config_petsc petscinc`
	     rm -fr config_petsc
	     have_petsc="yes"
	  fi		
	fi
	if test x"${have_petsc}" = xno; then
	   AC_MSG_RESULT(no)
	   PKG_CHECK_MODULES([pkgconfig_PETSc],[PETSc],
	   [have_petsc="yes"
	   PETSC_CPPFLAGS="$pkgconfig_PETSc_CFLAGS"
	   PETSC_LDFLAGS="$pkgconfig_PETSc_LIBS"
	   have_pkgconfig_petsc="yes"],
	   [have_petsc="no"])
	fi
	AC_SUBST(PETSC_LDFLAGS)
	if test x"${have_petsc}" = xyes; then
	   AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
	   CPPFLAGS="$CPPFLAGS $PETSC_CPPFLAGS"
	   LDFLAGS="$LDFLAGS $PETSC_LDFLAGS"
	   if test "x${ax_cv_cxx_compiler_vendor}" = xsgi; then	
	      LDFLAGS="$LDFLAGS -lfpe"
	   fi
	   if test x"${have_pkgconfig_petsc}" = xno; then
	      AC_MSG_RESULT(yes)
	   fi
	else
	   CPPFLAGS="$CPPFLAGS_SAVED"
	   LDFLAGS="$LDFLAGS_SAVED"
           AC_MSG_RESULT(no)
	fi
])



 
