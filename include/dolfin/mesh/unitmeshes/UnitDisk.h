// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UNIT_DISK_H
#define __DOLFIN_UNIT_DISK_H

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/// Triangular mesh of the 2D unit circle.
/// Given the number of cells (nx, ny) in each direction,
/// the total number of triangles will be 2*nx*ny and the
/// total number of vertices will be (nx + 1)*(ny + 1).

/// The Type is an enumerater taking values in {left, right or
/// crisscross} indicating the direction of the diagonals for
/// left/right or both == crisscross. The default is right.

class UnitDisk : public Mesh
{
public:
  enum Type
  {
    right,
    left,
    crisscross
  };
  enum Transformation
  {
    maxn,
    sumn,
    rotsumn
  };

  UnitDisk( size_t         nx,
            Type           type           = crisscross,
            Transformation transformation = rotsumn );

private:
  auto transformx( real x, real y, Transformation transformation ) -> real;

  auto transformy( real x, real y, Transformation transformation ) -> real;

  inline auto max( real x, real y ) -> real
  {
    return ( ( x > y ) ? x : y );
  }
};

}

#endif
