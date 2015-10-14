// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2015.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#ifndef __DOLFIN_INTERVAL_CELL_H
#define __DOLFIN_INTERVAL_CELL_H

#include <dolfin/mesh/CellType.h>

namespace dolfin
{

/**
 *  @class  IntervalCell
 *
 *  @brief  This class implements functionality for interval meshes.
 *
 */

class IntervalCell : public CellType
{
  // UFC: Topological Dimension
  static uint const TD = 1;

  // UFC: Number of Entities
  static uint const NE[2];

  // UFC: Number of Vertices (per entity)
  static uint const NV[2];

  // UFC: Vertex Coordinates
  static real const VC[2][1];

  // UFC: Edge - Incident Vertices
  static uint const EIV[1][2];

public:

  /// Specify cell type and facet type
  IntervalCell();

  ///
  ~IntervalCell();

  /// Clone pattern
  CellType* clone() const { return new IntervalCell(*this); }

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

  /// Order entities locally (connectivity 1-0)
  void orderEntities(Cell& cell) const;

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

  /// Compute (generalized) volume (length) of interval
  real volume(MeshEntity const& entity) const;

  /// Compute diameter of interval
  real diameter(MeshEntity const& entity) const;

  /// Compute circumradius of interval
  real circumradius(MeshEntity const& entity) const;

  /// Compute coordinates of midpoint
  Point midpoint(MeshEntity const& entity) const;

  /// Compute of given facet with respect to the cell
  Point normal(Cell const& cell, uint facet) const;

  /// Compute the area/length of given facet with respect to the cell
  real facetArea(Cell const& cell, uint facet) const;

  /// Check if point p intersects the entity
  bool intersects(MeshEntity const& e, Point const& p) const;

  /// Check if points line connecting p1 and p2 cuts the entity
  bool intersects(MeshEntity const& e, Point const& p1, Point const& p2) const;

  /// Create a mesh consisting of the reference cell
  Mesh create_reference_cell() const;

  /// Return description of cell type
  std::string description() const;

  /// Display information
  void disp() const;

  /// Check
  void check(Cell& cell) const;

};

} /* namespace dolfin */

#endif /* __DOLFIN_INTERVAL_CELL_H */
