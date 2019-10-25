// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

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

  /// Byteswap single precision floating point data
  float bswap(float x);

  /// Byteswap an array
  template<uint N, typename T> inline void bswap(T x[])
  {
    for (uint i = 0; i < N; i++) { x[i] = bswap(x[i]); }
  }

  /// Byteswap an array
  template<typename T> inline void bswap(T * begin, T * end)
  {
    while (begin != end) { *begin = bswap(*begin); ++begin; }
  }

}

#endif
