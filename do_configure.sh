#!/bin/sh

ARCH=amd64
if [ "X" = "X${CTLOPTROOT}" ];
then
	echo "CTLOPTROOT not set.";
	exit;
fi;

#export gts_CFLAGS="-I${CTLOPTROOT}/include -Wl,-rpath ${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib -R${CTLOPTROOT}/lib"
export gts_LIBS="/usr/lib/libgts.so"
#export gts_LIBS="-Wl,-rpath /usr/lib/libgts.so"
export CXX=mpicxx
export CXXFLAGS=" -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib "
export LDFLAGS=" -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib "

# --enable-optimize-p1 \

./configure --prefix=${CTLOPTROOT} \
	--libdir=${CTLOPTROOT}/lib \
	--enable-shared --with-pic \
	--enable-function-cache \
	--disable-boost-tr1 \
	--with-petsc --with-petsc-dir=${CTLOPTROOT} \
	--enable-mpi --enable-mpi-io \
	--enable-ufl \
	--with-parmetis-libdir=${CTLOPTROOT}/lib \
	--with-gts \
	--with-boost=/usr/include \
	--disable-progress-bar

