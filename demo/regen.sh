#!/bin/sh

for i in `find . -name *.ufl`; do pushd `dirname $i`;  ffc -l dolfin `basename $i`; popd; done;
