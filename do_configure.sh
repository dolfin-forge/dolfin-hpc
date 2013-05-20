#!/bin/sh

ARCH=amd64
if [ "X" == "X${CTLOPTROOT}" ];
then
	echo "CTLOPTROOT not set.";
	exit;
fi;

#export gts_CFLAGS="-I${CTLOPTROOT}/include -Wl,-rpath ${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib -R${CTLOPTROOT}/lib"
export gts_LIBS="${CTLOPTROOT}/lib/libgts.so"
#export gts_LIBS="-Wl,-rpath /usr/lib/libgts.so"
export CXX=mpicxx
export CXXFLAGS="-m64 -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib -R${CTLOPTROOT}/lib -L/usr/lib/amd64/atlas"
export LDFLAGS="-m64 -Wl,-rpath,${CTLOPTROOT}/lib -L${CTLOPTROOT}/lib -R${CTLOPTROOT}/lib -L/usr/lib/amd64/atlas"
export LIBS="-L/usr/lib/amd64/atlas -lcblas -lf77blas -llapack -latlas -lgfortran"

# --enable-optimize-p1 \

./configure --prefix=${CTLOPTROOT} \
	--libdir=${CTLOPTROOT}/lib \
	--enable-shared --with-pic \
	--enable-function-cache \
	--disable-boost-tr1 \
	--with-petsc --with-petsc-dir=${CTLOPTROOT} \
	--enable-mpi --enable-mpi-io \
	--enable-quadrature \
	--enable-ufl \
	--with-parmetis-libdir=${CTLOPTROOT}/lib \
	--with-gts \
	--with-boost=/usr/g++ \
	--disable-progress-bar

