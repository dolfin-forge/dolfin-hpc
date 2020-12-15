// Copyright (C) 2015 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_HEXAHEDRON_CELL_H
#define __DOLFIN_HEXAHEDRON_CELL_H

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshGeometry.h>

namespace dolfin
{

/**
 *  @class  HexahedronCell
 *
 *  @brief  This class implements functionality for hexahedral meshes.
 *
 */

class HexahedronCell : public CellType
{
  // UFC: Topological Dimension
  static uint const TD = 3;

  // UFC: Number of Entities
  static uint const NE[4][4];

  // UFC: Vertex Coordinates
  static real const VC[8][3];

  // UFC: Edge - Incident Vertices
  static uint const EIV[12][2];

  // UFC: Edge - Non-Incident Vertices
  static uint const ENV[12][6];

  // UFC: Face - Incident Vertices
  static uint const FIV[6][4];

  // UFC: Face - Non-Incident Vertices
  static uint const FNV[6][4];

public:

  /// Specify cell type and facet type
  HexahedronCell();

  ///
  ~HexahedronCell() override = default;

  /// Clone pattern
  CellType* clone() const override { return new HexahedronCell(*this); }

  /// Return topological dimension of cell
  uint dim() const override;

  /// Return number of entitites of given topological dimension
  uint num_entities(uint dim) const override;

  /// Return number of entities of given topological dimensions
  uint num_entities(uint d0, uint d1) const override;

  /// Return number of vertices for entity of given topological dimension
  uint num_vertices(uint dim) const override;

  /// Return orientation of the cell
  uint orientation(Cell const& cell) const override;

  /// Create entities e of given topological dimension from vertices v
  void create_entities(uint** e, uint dim, uint const* v) const override;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1, 3-0, 3-1, 3-2)
  void order_entities(MeshTopology& topology, uint i) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet(uint vertices[], Facet& facet) const override;

  /// Return if mesh connectivities require ordering
  bool connectivity_needs_ordering(uint d0, uint d1) const override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities(Mesh& mesh) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const override;

  /// Number of cells created by refinement pattern
  uint num_refined_cells() const override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  uint num_refined_vertices(uint dim) const override;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (area) of hexahedron
  real volume(MeshEntity const& entity) const override;

  /// Compute diameter of hexahedron
  real diameter(MeshEntity const& entity) const override;

  /// Compute circumradius of hexahedron
  real circumradius(MeshEntity const& entity) const override;

  /// Compute inradius of hexahedron
  real inradius(MeshEntity const& entity) const override;

  /// Compute coordinates of midpoint
  void midpoint(MeshEntity const& entity, real * p) const override;

  /// Compute of given facet with respect to the cell
  void normal(Cell const& cell, uint facet, real * n) const override;

  /// Compute the area/length of given facet with respect to the cell
  real facet_area(Cell const& cell, uint facet) const override;

  /// Check if point p intersects the entity
  bool intersects(MeshEntity const& e, Point const& p) const override;

  /// Check if points line connecting p1 and p2 cuts the entity
  bool intersects(MeshEntity const& e, Point const& p1, Point const& p2) const override;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  void create_reference_cell(Mesh& mesh) const override;

  /// Return coordinates of vertices in the reference cell
  real const * reference_vertex(uint i) const override;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  std::string description() const override;

  /// Display information
  void disp() const override;

  /// Check
  bool check(Cell& cell) const override;

private:

  // Find local index of edge i according to ordering convention
  uint findEdge(uint i, Cell const& cell) const;

  // Find local index of face i according to ordering convention
  uint findFace(uint i, Cell const& cell) const;

};

//-----------------------------------------------------------------------------
inline uint HexahedronCell::dim() const
{
  return 3;
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[3][dim];
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::num_entities(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD);
  dolfin_assert(d1 <= TD);
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get the coordinates of vertices v0, v1, v3 and v4
  MeshGeometry const& geometry = cell.mesh().geometry();
  Array<uint> const & vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v3 = geometry.x(vertices[3]);
  real const * v4 = geometry.x(vertices[4]);

  // Check whether (v0v1, v0v3, v0v4) is counter-clockwise
  real const S =
    + ((v1[1] - v0[1]) * (v3[2] - v0[2]) - (v1[2] - v0[2]) * (v3[1] - v0[1]))
          * (v4[0] - v0[0])
    + ((v1[2] - v0[2]) * (v3[0] - v0[0]) - (v1[0] - v0[0]) * (v3[2] - v0[2]))
          * (v4[1] - v0[1])
    + ((v1[0] - v0[0]) * (v3[1] - v0[1]) - (v1[1] - v0[1]) * (v3[0] - v0[0]))
          * (v4[2] - v0[2]);
  return (S < 0.0 ? 1 : 0);
}

//-----------------------------------------------------------------------------
inline bool HexahedronCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  // Do not order cell - vertices connectivities
  return (d0 > 0 && d0 > d1) && !(d0 == TD && d1 == 0);
}

//-----------------------------------------------------------------------------
inline void HexahedronCell::initialize_connectivities(Mesh& mesh) const
{
  mesh.init(1, 0);
  mesh.init(2, 0);
  mesh.init(2, 1);
  mesh.init(3, 0);
  mesh.init(3, 1);
  mesh.init(3, 2);
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::num_refined_cells() const
{
  return 8;
}

//-----------------------------------------------------------------------------
inline uint HexahedronCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  MAYBE_UNUSED(dim);
  return 1;
}

//-----------------------------------------------------------------------------
inline real HexahedronCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v3 = geometry.x(vertices[3]);
  real const * v4 = geometry.x(vertices[4]);

  // Compute mixed product of (v0v1, v0v3, v0v4)
  real const S =
      + ((v1[1] - v0[1]) * (v3[2] - v0[2]) - (v1[2] - v0[2]) * (v3[1] - v0[1]))
            * (v4[0] - v0[0])
      + ((v1[2] - v0[2]) * (v3[0] - v0[0]) - (v1[0] - v0[0]) * (v3[2] - v0[2]))
            * (v4[1] - v0[1])
      + ((v1[0] - v0[0]) * (v3[1] - v0[1]) - (v1[1] - v0[1]) * (v3[0] - v0[0]))
            * (v4[2] - v0[2]);
  return std::fabs(S);
}

//-----------------------------------------------------------------------------
inline real HexahedronCell::diameter(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);

  // Compute maximum diagonal
  real d0 = 0.0;
  real d1 = 0.0;
  real d2 = 0.0;
  real d3 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    d0 += (x6[i] - x0[i]) * (x6[i] - x0[i]);
    d1 += (x7[i] - x1[i]) * (x7[i] - x1[i]);
    d2 += (x4[i] - x2[i]) * (x4[i] - x2[i]);
    d3 += (x5[i] - x3[i]) * (x5[i] - x3[i]);
  }
  return std::sqrt(std::max(std::max(d0, d1),std::max(d2, d3)));
}

//-----------------------------------------------------------------------------
inline real HexahedronCell::circumradius(MeshEntity const& entity) const
{
  /// @todo No better idea for now
  return 0.5*std::sqrt(3.0)*this->diameter(entity);
}

//-----------------------------------------------------------------------------
inline real HexahedronCell::inradius(MeshEntity const& entity) const
{
  /// @todo No better idea for now
  return 0.5*this->diameter(entity);
}

//-----------------------------------------------------------------------------
inline void HexahedronCell::midpoint(MeshEntity const& entity, real * p) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    p[i] = 0.125
        * (x0[i] + x1[i] + x2[i] + x3[i] + x4[i] + x5[i] + x6[i] + x7[i]);
  }
}

//-----------------------------------------------------------------------------
inline void HexahedronCell::normal(Cell const& cell, uint facet, real * n) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();

  // first non-incident vertex
  real const * vf = geometry.x(c.entities(0)[FNV[facet][0]]);

  // vertices on the facet
  Array< uint > const & vertices = f.entities(0);
  real const * v0 = geometry.x(vertices[FIV[facet][0]]);
  real const * v1 = geometry.x(vertices[FIV[facet][1]]);
  real const * v3 = geometry.x(vertices[FIV[facet][3]]);

  // Vector normal to facet
  n[0] = + (v1[1] - v0[1]) * (v3[2] - v0[2])
         - (v1[2] - v0[2]) * (v3[1] - v0[1]);
  n[1] = + (v1[2] - v0[2]) * (v3[0] - v0[0])
         - (v1[0] - v0[0]) * (v3[2] - v0[2]);
  n[2] = + (v1[0] - v0[0]) * (v3[1] - v0[1])
         - (v1[1] - v0[1]) * (v3[0] - v0[0]);
  real const nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] /= nn;
  n[1] /= nn;
  n[2] /= nn;

  // Flip direction of normal so it points outward
  if (+ n[0]*(vf[0]-v0[0]) + n[1]*(vf[1]-v0[1]) + n[2]*(vf[2] - v0[2]) > 0.0)
  {
    n[0] *= -1.0;
    n[1] *= -1.0;
    n[2] *= -1.0;
  }
}

//-----------------------------------------------------------------------------
inline real HexahedronCell::facet_area(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);

  // Get the coordinates of the two vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  real const * p = geometry.x(f.entities(0)[FIV[facet][0]]);
  real const * q = geometry.x(f.entities(0)[FIV[facet][1]]);
  real const * r = geometry.x(f.entities(0)[FIV[facet][2]]);
  real const * s = geometry.x(f.entities(0)[FIV[facet][3]]);

  // Compute the area with diagonal formula
  real c0 = (r[1] - p[1]) * (s[2] - q[2]) - (r[2] - p[2]) * (s[1] - q[1]);
  real c1 = (r[2] - p[2]) * (s[0] - q[0]) - (r[0] - p[0]) * (s[2] - q[2]);
  real c2 = (r[0] - p[0]) * (s[1] - q[1]) - (r[1] - p[1]) * (s[0] - q[0]);
  return 0.5 * std::sqrt(c0 * c0 + c1 * c1 + c2 * c2);
}

//-----------------------------------------------------------------------------
inline real const * HexahedronCell::reference_vertex(uint i) const
{
  return &VC[i][0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif  /* __DOLFIN_HEXAHEDRON_CELL_H */
