// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BYTESWAP_H
#define __DOLFIN_BYTESWAP_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// Byteswap integer data
int bswap( int x );

/// Byteswap (unsigned) integer data
// uint bswap( uint x );

/// Byteswap double precision floating point data
real bswap( real x );

/// Byteswap single precision floating point data
float bswap( float x );

/// Byteswap an array
template < uint N, typename T >
void bswap( T x[] );

/// Byteswap an array
template < typename T >
void bswap( T * begin, T * end );

//-----------------------------------------------------------------------------
inline int bswap( int x )
{
  union
  {
    int           x;
    unsigned char b[4];
  } ein, eout;

  ein.x     = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];

  return eout.x;
}
//-----------------------------------------------------------------------------
// inline uint bswap( uint x )
// {
//   union
//   {
//     uint          x;
//     unsigned char b[4];
//   } ein, eout;

//   ein.x     = x;
//   eout.b[0] = ein.b[3];
//   eout.b[1] = ein.b[2];
//   eout.b[2] = ein.b[1];
//   eout.b[3] = ein.b[0];

//   return eout.x;
// }
//-----------------------------------------------------------------------------
inline uint bswap( uint32_t x )
{
  union
  {
    unsigned int  x;
    unsigned char b[4];
  } ein, eout;

  ein.x     = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];

  return eout.x;
}
//-----------------------------------------------------------------------------
inline real bswap( real x )
{
  union
  {
    real          x;
    unsigned char b[8];
  } ein, eout;

  ein.x     = x;
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
inline float bswap( float x )
{
  union
  {
    float         x;
    unsigned char b[4];
  } ein, eout;

  ein.x     = x;
  eout.b[0] = ein.b[3];
  eout.b[1] = ein.b[2];
  eout.b[2] = ein.b[1];
  eout.b[3] = ein.b[0];
  return eout.x;
}
//-----------------------------------------------------------------------------
template < uint N, typename T >
inline void bswap( T x[] )
{
  for ( uint i = 0; i < N; i++ )
  {
    x[i] = bswap( x[i] );
  }
}
//-----------------------------------------------------------------------------
template < typename T >
inline void bswap( T * begin, T * end )
{
  while ( begin != end )
  {
    *begin = bswap( *begin );
    ++begin;
  }
}
//-----------------------------------------------------------------------------

}

#endif
