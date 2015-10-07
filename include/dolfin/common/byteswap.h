// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-07-25
// Last changed: 2015-07-25

#ifndef __DOLFIN_BYTESWAP_H
#define __DOLFIN_BYTESWAP_H

#include <dolfin/common/types.h>

namespace dolfin
{

  /// Byteswap integer data
  int bswap(int x);

  /// Byteswap (unsigned) integer data
  uint bswap(uint x);

  /// Byteswap double precision floating point data
  real bswap(real x);

  /// Byteswap an array
  template<typename T> void bswap(T x[], uint n);

}

#endif
