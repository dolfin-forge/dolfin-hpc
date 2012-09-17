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
	have_petsc=no
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
	     ac_petsc_libdir="$ac_petsc_dir/lib"
	     PETSC_PKG_LIBS=`grep PACKAGES_LIBS $ac_petsc_dir/conf/petscvariables | sed 's/PACKAGES_LIBS =/''/'`
	     PETSC_MPI_INC=`grep 'MPI_INCLUDE =' $ac_petsc_dir/conf/petscvariables | sed 's/MPI_INCLUDE =/''/'`
	     PETSC_CPPFLAGS="-I$ac_petsc_dir/include/ $PETSC_MPI_INC"
	     PETSC_SYS_LIB=`grep 'PETSC_SYS_LIB_BASIC ' ${ac_petsc_dir}/conf/variables | sed 's/PETSC_SYS_LIB_BASIC/''/' | sed 's/=/''/'`
	     PETSC_SYS_LIB=`echo $PETSC_SYS_LIB | sed 's/ //g'`

	     PETSC_VEC_LIB=`grep 'PETSC_VEC_LIB_BASIC ' $ac_petsc_dir/conf/variables | sed 's/PETSC_VEC_LIB_BASIC/''/' | sed 's/=/''/' | sed 's/${PETSC_SYS_LIB_BASIC}/''/'`
	     PETSC_VEC_LIB=`echo $PETSC_VEC_LIB | sed 's/ //g'`

	     PETSC_MAT_LIB=`grep 'PETSC_MAT_LIB_BASIC '  $ac_petsc_dir/conf/variables | sed 's/PETSC_MAT_LIB_BASIC/''/'| sed 's/=/''/' | sed 's/${PETSC_VEC_LIB_BASIC}/''/'`
	     PETSC_MAT_LIB=`echo $PETSC_MAT_LIB | sed 's/ //g'`

	     PETSC_DM_LIB=`grep 'PETSC_DM_LIB_BASIC ' $ac_petsc_dir/conf/variables | sed 's/PETSC_DM_LIB_BASIC/''/'| sed 's/=/''/' | sed 's/${PETSC_MAT_LIB_BASIC}/''/'`
	     PETSC_DM_LIB=`echo $PETSC_DM_LIB | sed 's/ //g'`

	     PETSC_KSP_LIB=`grep 'PETSC_KSP_LIB_BASIC ' $ac_petsc_dir/conf/variables | sed 's/PETSC_KSP_LIB_BASIC/''/' | sed 's/=/''/' | sed 's/${PETSC_DM_LIB_BASIC}/''/'`
	     PETSC_KSP_LIB=`echo $PETSC_KSP_LIB | sed 's/ //g'`

	     PETSC_SNES_LIB=`grep 'PETSC_SNES_LIB_BASIC ' $ac_petsc_dir/conf/variables | sed 's/PETSC_SNES_LIB_BASIC/''/' | sed 's/=/''/' |  sed 's/^[ \t]*//' | sed 's/${PETSC_KSP_LIB_BASIC}/''/'`
	     PETSC_SNES_LIB=`echo $PETSC_SNES_LIB | sed 's/ //g'`

	     PETSC_TS_LIB=`grep 'PETSC_TS_LIB_BASIC ' $ac_petsc_dir/conf/variables | sed 's/PETSC_TS_LIB_BASIC/''/' | sed 's/=/''/' | sed 's/${PETSC_SNES_LIB_BASIC}/''/'`
	     PETSC_TS_LIB=`echo $PETSC_TS_LIB | sed 's/ //g'`

	     PETSC_LDFLAGS="-L$ac_petsc_libdir $PETSC_PKG_LIBS $PETSC_TS_LIB $PETSC_SNES_LIB $PETSC_KSP_LIB $PETSC_DM_LIB $PETSC_MAT_LIB $PETSC_VEC_LIB $PETSC_SYS_LIB"
#	     PETSC_LDFLAGS=`echo $PETSC_LDFLAGS | sed 's/ //g'`
	     have_petsc="yes"
	   fi
		
	fi
	AC_SUBST(PETSC_LDFLAGS)
	if test x"${have_petsc}" = xyes; then
	   AC_DEFINE(HAVE_PETSC,1,[Define if you have the Petsc library.])
	   CPPFLAGS="$CPPFLAGS $PETSC_CPPFLAGS"
	   LDFLAGS="$LDFLAGS $PETSC_LDFLAGS"
	   if test "x${ax_cv_cxx_compiler_vendor}" = xsgi; then	
	      LDFLAGS="$LDFLAGS -lfpe"
	   fi
	   AC_MSG_RESULT(yes)
	else
	   CPPFLAGS="$CPPFLAGS_SAVED"
	   LDFLAGS="$LDFLAGS_SAVED"
           AC_MSG_RESULT(no)
	fi
])



 