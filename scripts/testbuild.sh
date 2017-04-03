#!/bin/sh

# This script is just a convenient way to run test builds in a streamlined
# fashion for Continuous Integration e.g. with Jenkins.

if [[ "X${COMPILER}" = "Xgcc" ]]
then
  CXXFLAGS="${CXXFLAGS} -O3 -pedantic -Wno-long-long -Werror=reorder"
  CXXFLASG="${CXXFLAGS} -Wno-unused-function -Werror=return-local-addr"
  CXXFLASG="${CXXFLAGS} -fno-fast-math -finline -fbounds-check"
fi

#-----------------------------------------------------------------------------
config() {
  if [[ "X$1" = "Xtrue" ]]
  then CONFIGURE_OPTIONS="${CONFIGURE_OPTIONS} $2"
  fi
}

CONFIGURE_OPTIONS="--prefix=${PWD}/build \
	           --enable-static --disable-shared --with-pic \
                   --enable-check --enable-debug"

config "${MPI}"		"--enable-mpi"
echo "MPI       : ${MPI}"
config "${MPI_IO}"	"--enable-mpi-io"
echo "MPI_IO    : ${MPI_IO}"
config "${OPENMP}"	"--enable-openmp"
echo "OPENMP    : ${OPENMP}"
config "${GTS}"		"--with-gts"
echo "GTS       : ${GTS}"
config "${JANPACK}"	"--with-janpack"
echo "JANPACK   : ${JANPACK}"
config "${PETSC}"	"--with-petsc=${PETSC_DIR}"
echo "PETSC     : ${PETSC}"
config "${PARMETIS}"	"--with-parmetis=${PARMETIS_DIR}"
echo "PARMETIS  : ${PARMETIS}"
config "${XML}"		"--with-xml"
echo "XML       : ${XML}"

#-----------------------------------------------------------------------------
NJOBS=${NJOBS:-1}
echo "Build on ${NJOBS} jobs"

fresh_build() {
  ( sh regen.sh && \
    ./configure ${CONFIGURE_OPTIONS} CXXFLAGS="$CXXFLAGS" && \
    make -j "${NJOBS}" clean check dist )
}

fresh_build && \
tar xvzf dolfin-0.9.0-hpc.tar.gz && \
sh regen.sh && \
fresh_build

#-----------------------------------------------------------------------------
