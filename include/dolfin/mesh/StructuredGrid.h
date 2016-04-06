// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_STRUCTURED_GRID_H
#define __DOLFIN_STRUCTURED_GRID_H

#include <dolfin/mesh/Mesh.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/BoundingBox.h>
#include <dolfin/mesh/CellType.h>

namespace dolfin
{

/**
 *  @class  StructuredGrid
 *
 *  @brief  Implements the generation of structured grids
 *
 */

class StructuredGrid : public Mesh
{

public:

  /// Create unit structured grid
  StructuredGrid(CellType const& type, uint N);

  /// Destructor
  ~StructuredGrid();

private:

  //
  void init(CellType const& type);

  //--- ATTRIBUTES ------------------------------------------------------------

  BoundingBox bbox_;
  uint const n_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_STRUCTURED_GRID_H */
