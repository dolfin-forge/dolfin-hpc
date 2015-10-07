// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
// Modified by Garth N. Wells 2006.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#ifndef __DOLFIN_TETRAHEDRON_CELL_H
#define __DOLFIN_TETRAHEDRON_CELL_H

#include "CellType.h"

namespace dolfin
{

class Cell;

/// This class implements functionality for tetrahedral meshes.

class TetrahedronCell : public CellType
{
public:

  /// Specify cell type and facet type
  TetrahedronCell();

  /// Return topological dimension of cell
  uint dim() const;

  /// Return number of entitites of given topological dimension
  uint numEntities(uint dim) const;

  /// Return number of vertices for entity of given topological dimension
  uint numVertices(uint dim) const;

  /// Return orientation of the cell
  uint orientation(Cell const& cell) const;

  /// Create entities e of given topological dimension from vertices v
  void createEntities(uint** e, uint dim, uint const* v) const;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1, 3-0, 3-1, 3-2)
  void orderEntities(Cell& cell) const;

  /// Regular refinement of cell
  void refineCell(Cell& cell, MeshEditor& editor, uint& current_cell) const;

  /// Irregular refinement of cell
  void refineCellIrregular(Cell& cell, MeshEditor& editor, uint& current_cell,
                           uint refinement_rule, uint* marked_edges) const;

  /// Compute volume of tetrahedron
  real volume(MeshEntity const& tetrahedron) const;

  /// Compute diameter of tetrahedron
  real diameter(MeshEntity const& tetrahedron) const;

  /// Compute circumradius of tetrahedron
  real circumradius(MeshEntity const& tetrahedron) const;

  /// Compute component i of normal of given facet with respect to the cell
  real normal(Cell const& cell, uint facet, uint i) const;

  /// Compute of given facet with respect to the cell
  Point normal(Cell const& cell, uint facet) const;

  /// Compute the area/length of given facet with respect to the cell
  real facetArea(Cell const& cell, uint facet) const;

  /// Check if point p intersects the cell
  bool intersects(MeshEntity const& entity, Point const& p) const;

  /// Check if points line connecting p1 and p2 cuts the cell
  bool intersects(MeshEntity const& entity, Point const& p1,
                  Point const& p2) const;

  /// Return description of cell type
  std::string description() const;

private:

  // Find local index of edge i according to ordering convention
  uint findEdge(uint i, Cell const& cell) const;

};

}

#endif
