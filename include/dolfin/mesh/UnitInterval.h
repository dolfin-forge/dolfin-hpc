// Copyright (C) 2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-23
// Last changed: 2007-11-23

#ifndef __DOLFIN_UNIT_INTERVAL_H
#define __DOLFIN_UNIT_INTERVAL_H

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/// Interval mesh of the 1D unit line (0,1).
/// Given the number of cells (nx) in the axial direction,
/// the total number of intervals will be nx and the
/// total number of vertices will be (nx + 1).

class UnitInterval : public Mesh
{

public:

  UnitInterval(uint nx);

};

} /* namespace dolfin */

#endif /* __DOLFIN_UNIT_INTERVAL_H */
