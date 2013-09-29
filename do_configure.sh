#!/bin/sh

ARCH=amd64
if [ "X" = "X${CTLOPTROOT}" ];
then
	echo "CTLOPTROOT not set.";
fi;

#export gts_LIBS="/usr/lib/libgts.so"
export CXX=mpicxx
export CFLAGS="-m64 -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib "
export CXXFLAGS="-m64 -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib -g -std=c++0x -Wall -pedantic-errors -fno-fast-math -finline -fbounds-check "
export LDFLAGS="-m64 -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib "

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

