// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Aurelien Larcher, 2015.
//
// First added:  2007-12-12
// Last changed: 2008-06-20

#ifndef __DOLFIN_POINT_CELL_H
#define __DOLFIN_POINT_CELL_H

#include <dolfin/mesh/CellType.h>

namespace dolfin
{

/**
 *  @class  PointCell
 *
 *  @brief  This class implements functionality for point meshes.
 *
 */

class PointCell : public CellType
{
  // UFC: Topological Dimension
  static uint const TD = 0;

  // UFC: Number of Entities
  static uint const NE[1];

  // UFC: Number of Vertices (per entity)
  static uint const NV[1];

  // UFC: Vertex Coordinates
  static real const VC[1][1];

public:

  /// Specify cell type and facet type
  PointCell();

  ///
  ~PointCell();

  /// Clone pattern
  CellType* clone() const { return new PointCell(*this); }

  /// Return topological dimension of cell
  uint dim() const;

  /// Return number of entitites of given topological dimension
  uint num_entities(uint dim) const;

  /// Return number of vertices for entity of given topological dimension
  uint num_vertices(uint dim) const;

  /// Return orientation of the cell
  uint orientation(Cell const& cell) const;

  /// Create entities e of given topological dimension from vertices v
  void create_entities(uint** e, uint dim, uint const* v) const;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1)
  void order_entities(Cell& cell) const;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const;

  /// Number of cells created by refinement pattern
  uint num_refined_cells() const;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  uint num_refined_vertices(uint dim) const;

  /// Return if refinement pattern requires entities of given dimension
  bool refinement_needs_entities(uint dim) const;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (area) of triangle
  real volume(MeshEntity const& entity) const;

  /// Compute diameter of triangle
  real diameter(MeshEntity const& entity) const;

  /// Compute circumradius of triangle
  real circumradius(MeshEntity const& entity) const;

  /// Compute coordinates of midpoint
  Point midpoint(MeshEntity const& entity) const;

  /// Compute of given facet with respect to the cell
  Point normal(Cell const& cell, uint facet) const;

  /// Compute the area/length of given facet with respect to the cell
  real facet_area(Cell const& cell, uint facet) const;

  /// Check if point p intersects the cell
  bool intersects(MeshEntity const& e, Point const& p) const;

  /// Check if points line connecting p1 and p2 cuts the cell
  bool intersects(MeshEntity const& e, Point const& p1, Point const& p2) const;

  /// Create a mesh consisting of the reference cell
  Mesh create_reference_cell() const;

  /// Return description of cell type
  std::string description() const;

  /// Display information
  void disp() const;

  /// Check
  void check(Cell& cell) const;

private:

};

} /* namespace dolfin */

#endif /* __DOLFIN_POINT_CELL_H */
