#!/bin/sh


echo "Updating configuration..."
echo "Running libtoolize"
if which libtoolize > /dev/null 2>&1; then
    libtoolize -i
elif which glibtoolize > /dev/null 2>&1; then
    glibtoolize -i
else
    echo "No libtoolize found on your system"
    exit 1
fi

echo "Running aclocal"
if which aclocal > /dev/null 2>&1; then
if test -n "${NO_ACLOCAL_INSTALL}"; then
  aclocal -I m4
else
  aclocal -I m4 --install
fi
else
  echo "No aclocal found on your system"
  exit 1
fi

echo "Running autoconf"
autoconf
echo "Running automake"
automake -a

echo "Deleting autom4te.cache directory"
rm -r autom4te.cache

