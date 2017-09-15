// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-07-25
// Last changed: 2015-07-25

#include <dolfin/common/byteswap.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
int bswap(int x)
{
  union {
    int x;
    unsigned char b[4];
  } ein, eout;

  ein.x = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];

  return eout.x;
}
//-----------------------------------------------------------------------------
uint bswap(uint x)
{
  union {
    uint x;
    unsigned char b[4];
  } ein, eout;

  ein.x = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];

  return eout.x;
}
//-----------------------------------------------------------------------------
real bswap(real x)
{
  union {
    real x;
    unsigned char b[8];
  } ein, eout;

  ein.x = x;
  eout.b[0] = ein.b[7];
  eout.b[1] = ein.b[6];
  eout.b[2] = ein.b[5];
  eout.b[3] = ein.b[4];
  eout.b[4] = ein.b[3];
  eout.b[5] = ein.b[2];
  eout.b[6] = ein.b[1];
  eout.b[7] = ein.b[0];
  return eout.x;
}
//-----------------------------------------------------------------------------
float bswap(float x)
{
  union {
    float x;
    unsigned char b[4];
  } ein, eout;

  ein.x = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];
  return eout.x;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
