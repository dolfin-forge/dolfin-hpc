// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
// Modified by Kristian Oelgaard, 2006-2007.
// Modified by Dag Lindbo, 2008.
// Modified by Aurelien Larcher, 2014-2015.
//
// First added:  2006-06-05
// Last changed: 2014-11-07

#include <dolfin/mesh/TriangleCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/GeometricPredicates.h>
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const TriangleCell::NE[3] =
{ 3, 3, 1 };

// UFC: Number of Vertices (per entity)
uint const TriangleCell::NV[3] =
{ 1, 2, 3 };

// UFC: Vertex Coordinates
real const TriangleCell::VC[3][2] =
{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };

// UFC: Edge - Incident Vertices
uint const TriangleCell::EIV[3][2] =
{ { 1, 2 }, { 0, 2 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
uint const TriangleCell::ENV[3][1] =
{ { 0 }, { 1 }, { 2 } };

//-----------------------------------------------------------------------------
TriangleCell::TriangleCell() :
    CellType("triangle", CellType::triangle, CellType::interval)
{
}
//-----------------------------------------------------------------------------
TriangleCell::~TriangleCell()
{
}
//-----------------------------------------------------------------------------
uint TriangleCell::dim() const
{
  return 2;
}
//-----------------------------------------------------------------------------
uint TriangleCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim];
}
//-----------------------------------------------------------------------------
uint TriangleCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NV[dim];
}
//-----------------------------------------------------------------------------
uint TriangleCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v2 = geometry.x(vertices[2]);
  return (
      ((v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0]))
          < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void TriangleCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges
  if (dim != 1)
  {
    error("Invalid topological dimension for creation of entities: %d.", dim);
  }

  // Create the three edges following UFC convention
  e[0][0] = v[1];
  e[0][1] = v[2];
  e[1][0] = v[0];
  e[1][1] = v[2];
  e[2][0] = v[0];
  e[2][1] = v[1];
}
//-----------------------------------------------------------------------------
void TriangleCell::order_entities(MeshTopology& topology, uint i) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1
  dolfin_assert(topology.type(i).cellType() == this->cell_type);

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology.is_computed(1, 0))
  {
    dolfin_assert(topology.is_computed(2, 1));

    // Get edges
    uint* cell_edges = topology(2, 1)(i);

    // Sort vertices on each edge
    for (uint i = 0; i < 3; ++i)
    {
      uint* edge_vertices = topology(1, 0)(cell_edges[i]);
      std::sort(edge_vertices, edge_vertices + 2);
    }
  }

  // Sort local vertices on cell in ascending order, connectivity 2 - 0
  if (topology.is_computed(2, 0))
  {
    uint* cell_vertices = topology(2, 0)(i);
    std::sort(cell_vertices, cell_vertices + 3);
  }

  // Sort local edges on cell after non-incident vertex, connectivity 2 - 1
  if (topology.is_computed(2, 1))
  {
    dolfin_assert(topology.is_computed(2, 1));

    // Get cell vertices and edges
    uint* cell_vertices = topology(2, 0)(i);
    uint* cell_edges = topology(2, 1)(i);

    // Loop over vertices on cell
    for (uint i = 0; i < 3; ++i)
    {
      // Loop over edges on cell
      for (uint j = i; j < 3; ++j)
      {
        uint* edge_vertices = topology(1, 0)(cell_edges[j]);

        // Check if the ith vertex of the cell is non-incident with edge j
#if __SUNPRO_CC
        int n1 = 0;
        std::count(edge_vertices, edge_vertices + 2, cell_vertices[i], n1);
        if ( n1 == 0)
#else
        if (std::count(edge_vertices, edge_vertices + 2, cell_vertices[i]) == 0)
#endif
        {
          // Swap edge numbers
          uint tmp = cell_edges[i];
          cell_edges[i] = cell_edges[j];
          cell_edges[j] = tmp;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void TriangleCell::order_facet(uint vertices[], Facet& facet) const
{
  // Get mesh
  Mesh& mesh = facet.mesh();

  // Get the vertex opposite to the facet (the one we remove)
  uint vertex = 0;
  const Cell cell(mesh, facet.entities(mesh.topology_dimension())[0]);
  for (uint i = 0; i < cell.num_entities(0); i++)
  {
    bool not_in_facet = true;
    vertex = cell.entities(0)[i];
    for (uint j = 0; j < facet.num_entities(0); j++)
    {
      if (vertex == facet.entities(0)[j])
      {
        not_in_facet = false;
        break;
      }
    }
    if (not_in_facet)
    {
      break;
    }
  }

  // Order
  Point const p = mesh.geometry().point(vertex);
  Point p0 = mesh.geometry().point(facet.entities(0)[0]);
  Point p1 = mesh.geometry().point(facet.entities(0)[1]);
  Point v = p1 - p0;
  Point n(v[1], -v[0]);
  if (n.dot(p0 - p) < 0.0)
  {
    std::swap(vertices[0], vertices[1]);
  }
}
//-----------------------------------------------------------------------------
bool TriangleCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  return (d0 > 0 && d0 > d1);
}
//-----------------------------------------------------------------------------
void TriangleCell::initialize_connectivities(Mesh& mesh) const
{
  mesh.init(1, 0);
  mesh.init(2, 0);
  mesh.init(2, 1);
}
//-----------------------------------------------------------------------------
void TriangleCell::refine_cell(Cell& cell, MeshEditor& editor,
                               uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  uint const * v = cell.entities(0);
  dolfin_assert(v);
  uint const * e = cell.entities(1);
  dolfin_assert(e);

  // Compute indices for the six new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const offset = cell.mesh().size(0);
  uint const e0 = offset + e[findEdge(0, cell)];
  uint const e1 = offset + e[findEdge(1, cell)];
  uint const e2 = offset + e[findEdge(2, cell)];

  // Add the four new cells
  uint cv0[3] = { v0, e2, e1 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint cv1[3] = { v1, e0, e2 };
  editor.add_cell(current_cell++, &cv1[0]);
  uint cv2[3] = { v2, e1, e0 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint cv3[3] = { e0, e1, e2 };
  editor.add_cell(current_cell++, &cv3[0]);
}
//-----------------------------------------------------------------------------
uint TriangleCell::num_refined_cells() const
{
  return 4;
}
//-----------------------------------------------------------------------------
uint TriangleCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return (dim > 1 ? 0 : 1);
}
//-----------------------------------------------------------------------------
real TriangleCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);

  switch (geometry.dim())
    {
    case 2:
      // Compute area of triangle embedded in R^2
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5
          * std::abs(
              (x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1])
                  - (x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1]));
      break;
    case 3:
      // Compute area of triangle embedded in R^3
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5 * std::sqrt(
      + std::pow( (x0[1] * x1[2] + x0[2] * x2[1] + x1[1] * x2[2])
                    - (x2[1] * x1[2] + x2[2] * x0[1] + x1[1] * x0[2]), 2)
      + std::pow( (x0[2] * x1[0] + x0[0] * x2[2] + x1[2] * x2[0])
                    - (x2[2] * x1[0] + x2[0] * x0[2] + x1[2] * x0[0]), 2)
      + std::pow( (x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1])
                    - (x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1]), 2)
      );
      break;
    default:
      error("Volume of triangle only implemented for R^2 or R^3.");
      break;
    }
  return 0.0;
}
//-----------------------------------------------------------------------------
real TriangleCell::diameter(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    e1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    e2 += (x0[i] - x2[i]) * (x0[i] - x2[i]);
  }
  return std::sqrt(std::max(std::max(e0, e1), e2));
}
//-----------------------------------------------------------------------------
real TriangleCell::circumradius(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    e1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    e2 += (x0[i] - x2[i]) * (x0[i] - x2[i]);
  }
  e0 = std::sqrt(e0);
  e1 = std::sqrt(e1);
  e2 = std::sqrt(e2);

  //
  return 0.5
      * std::sqrt(
          (e0 + e1 - e2) * (e2 + e0 - e1) * (e1 + e2 - e0) / (e0 + e1 + e2));
}
//-----------------------------------------------------------------------------
real TriangleCell::inradius(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    e1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    e2 += (x0[i] - x2[i]) * (x0[i] - x2[i]);
  }
  e0 = std::sqrt(e0);
  e1 = std::sqrt(e1);
  e2 = std::sqrt(e2);

  // Formula for circumradius from http://mathworld.wolfram.com
  // Using Heron's formula for the volume instead of calling volume()
  real const s = 0.5 * (e0 + e1 + e2);
  return 0.25 * e0 * e1 * e2 / std::sqrt(s * (s - e0) * (s - e1) * (s - e2));
}
//-----------------------------------------------------------------------------
void TriangleCell::midpoint(MeshEntity const& entity, real * p) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[0]);

  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  uint const gdim = geometry.dim();
  for (uint i = 0; i < gdim; ++i)
  {
    p[i] = ( x0[i] + x1[i] + x2[i] ) / 3.0;
  }
}
//-----------------------------------------------------------------------------
void TriangleCell::normal(Cell const& cell, uint facet, real * n) const
{
  dolfin_assert(cell.type() == this->cell_type);

  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();
  // Get coordinates of opposite vertex
  real const * p0 = geometry.x(c.entities(0)[facet]);
  // Get coordinates of edge vertices
  uint const * vertices = f.entities(0);
  real const * p1 = geometry.x(vertices[0]);
  real const * p2 = geometry.x(vertices[1]);
  uint const gdim = geometry.dim();
  switch (gdim)
  {
    case 2:
      {
        n[0] = p2[1] - p1[1];
        n[1] = p1[0] - p2[0];
        real const nn = std::sqrt(n[0] * n[0] + n[1] * n[1]);
        n[0] /= nn;
        n[1] /= nn;
        // Flip direction of normal so it points outward
        if ((n[0] * (p1[0] - p0[0]) + n[1] * (p1[1] - p0[1])) < 0.0)
        {
          n[0] *= -1.0;
          n[1] *= -1.0;
        }
      }
      break;
    case 3:
      {
        real ac = (p2[0] - p1[0]) * (p2[0] - p0[0])
            + (p2[1] - p1[1]) * (p2[1] - p0[1])
            + (p2[2] - p1[2]) * (p2[2] - p0[2]);
        real ab = (p2[0] - p1[0]) * (p1[0] - p0[0])
            + (p2[1] - p1[1]) * (p1[1] - p0[1])
            + (p2[2] - p1[2]) * (p1[2] - p0[2]);
        n[0] = (p1[0] - p0[0]) * ac - (p2[0] - p0[0]) * ab;
        n[1] = (p1[1] - p0[1]) * ac - (p2[1] - p0[1]) * ab;
        n[2] = (p1[2] - p0[2]) * ac - (p2[2] - p0[2]) * ab;
        real nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        n[0] /= nn;
        n[1] /= nn;
        n[2] /= nn;
        dolfin_assert((n[0] * (p1[0] - p0[0]) + n[1] * (p1[1] - p0[1])
                        + n[2] * (p1[2] - p0[2])) > 0.0);
      }
      break;
    default:
      error("TriangleCell : facet normal not implemented in R^%u", gdim);
      break;
  }
}
//-----------------------------------------------------------------------------
real TriangleCell::facet_area(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = f.entities(0);
  real const * p0 = geometry.x(vertices[0]);
  real const * p1 = geometry.x(vertices[1]);
  // Compute distance between vertices
  real meas = 0.0;
  uint const gdim = geometry.dim();
  for (uint d = 0; d < gdim; ++d)
  {
    meas += (p1[d] - p0[d]) * (p1[d] - p0[d]);
  }
  return std::sqrt(meas);
}
//-----------------------------------------------------------------------------
bool TriangleCell::intersects(MeshEntity const& e, Point const& p) const
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[0]);

  // Get mesh geometry
  MeshGeometry const& geometry = e.mesh().geometry();

  // Get global index of vertices of the triangle
  uint const ort = orientation((Cell&) e);
  uint const v0 = e.entities(0)[0];
  uint const v1 = ( ort == 0 ? e.entities(0)[1] : e.entities(0)[2] );
  uint const v2 = ( ort == 0 ? e.entities(0)[2] : e.entities(0)[1] );

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x(v0);
  real const * x1 = geometry.x(v1);
  real const * x2 = geometry.x(v2);

  // The return value of orient2d is defined as twice the signed measure of
  // the triangle defined by the points given as arguments.
  real tol = geometry.abs_tolerance(2);
  real l0 = std::sqrt((x1[0] - x0[0])*(x1[0] - x0[0])
                        + (x1[1] - x0[1])*(x1[1] - x0[1]));
  real l1 = std::sqrt((x2[0] - x1[0])*(x2[0] - x1[0])
                        + (x2[1] - x1[1])*(x2[1] - x1[1]));
  real l2 = std::sqrt((x0[0] - x2[0])*(x0[0] - x2[0])
                        + (x0[1] - x2[1])*(x0[1] - x2[1]));

  if (geometry.dim() == 2)
  {
    // Test orientation of p w.r.t. each edge
    real d1 = orient2d( x0, x1, &p[0]);
    if (d1 < (- 2.0  * tol * l0)) return false;
    real d2 = orient2d( x1, x2, &p[0]);
    if (d2 < (- 2.0  * tol * l1)) return false;
    real d3 = orient2d( x2, x0, &p[0]);
    if (d3 < (- 2.0  * tol * l2)) return false;
  }
  else if  (geometry.dim() == 3)
  {
    real n0 = (x1[1] - x0[1])*(x2[2] - x0[2]) - (x1[2] - x0[2])*(x2[1] - x0[1]);
    real n1 = (x1[2] - x0[2])*(x2[0] - x0[0]) - (x1[0] - x0[0])*(x2[2] - x0[2]);
    real n2 = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
    real zz = (n2*n2) / (n0*n0 + n1*n1 + n2*n2);

    // Check direction of the normal to the triangle; |n.ez|/||n|| >= sqrt(2)/2
    if(zz >= 0.5)
    {
      // Test orientation of p w.r.t. each edge in (x,y) plane
      real d1 = orient2d( &x0[0], &x1[0], &p[0]);
      if (d1 < (- 2.0  * tol * l0)) return false;
      real d2 = orient2d( &x1[0], &x2[0], &p[0]);
      if (d2 < (- 2.0  * tol * l1)) return false;
      real d3 = orient2d( &x2[0], &x0[0], &p[0]);
      if (d3 < (- 2.0  * tol * l2)) return false;
    }
    else
    {
      // Test orientation of p w.r.t. each edge in (y,z) plane
      real d1 = orient2d( &x0[1], &x1[1], &p[1]);
      if (d1 < (- 2.0  * tol * l0)) return false;
      real d2 = orient2d( &x1[1], &x2[1], &p[1]);
      if (d2 < (- 2.0  * tol * l1)) return false;
      real d3 = orient2d( &x2[1], &x0[1], &p[1]);
      if (d3 < (- 2.0  * tol * l2)) return false;
    }
  }
  else
  {
    error("Collision of triangle with point only implemented in R^2 and R^3");
  }

  return true;
}
//-----------------------------------------------------------------------------
bool TriangleCell::intersects(MeshEntity const& e, Point const& p1,
                              Point const& p2) const
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[0]);

  // Get mesh geometry
  MeshGeometry const& geometry = e.mesh().geometry();

  // Get global index of vertices of the triangle
  uint const ort = orientation((Cell&) e);
  uint const v0 = e.entities(0)[0];
  uint const v1 = ( ort == 1 ? e.entities(0)[2] : e.entities(0)[1] );
  uint const v2 = ( ort == 1 ? e.entities(0)[1] : e.entities(0)[2] );

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x(v0);
  real const * x1 = geometry.x(v1);
  real const * x2 = geometry.x(v2);

  real d1, d2, d3;

  if (geometry.dim() == 2)
  {
    // Test orientation of each vertex w.r.t. p1-p2
    d1 = orient2d(&p1[0], &p2[0], x0);
    d2 = orient2d(&p1[0], &p2[0], x1);
    d3 = orient2d(&p1[0], &p2[0], x2);

    if (d1 < 0 && d2 < 0 && d3 < 0) return false;
    if (d1 > 0 && d2 > 0 && d3 > 0) return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x0-x1:
    d1 = orient2d(x0, x1, &p1[0]);
    d2 = orient2d(x0, x1, &p2[0]);

    if (d1 < 0 && d2 < 0) return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x1-x2:
    d1 = orient2d(x1, x2, &p1[0]);
    d2 = orient2d(x1, x2, &p2[0]);

    if (d1 < 0 && d2 < 0) return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x2-x0:
    d1 = orient2d(x2, x0, &p1[0]);
    d2 = orient2d(x2, x0, &p2[0]);

    if (d1 < 0 && d2 < 0) return false;
  }
  else
  {
    error("Collision of triangle with segment only implemented in R^2");
  }

  return true;
}
//-----------------------------------------------------------------------------
void TriangleCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::triangle, 2, DOLFIN_COMM_SELF);
  me.init_vertices(3);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.add_vertex(2, VC[2]);
  me.init_cells(1);
  uint const cv0[3] = { 0, 1, 2 };
  me.add_cell(0, cv0);
  me.close();
}
//-----------------------------------------------------------------------------
real const * TriangleCell::reference_vertex(uint i) const
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------
std::string TriangleCell::description() const
{
  return std::string("triangle (simplex of topological dimension 2)");
}
//-----------------------------------------------------------------------------
void TriangleCell::disp() const
{
  message("TriangleCell");
  begin(  "------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
uint TriangleCell::findEdge(uint i, Cell const& cell) const
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  uint const* v = cell.entities(0);
  dolfin_assert(v);
  uint const* e = cell.entities(1);
  dolfin_assert(e);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  for (uint j = 0; j < 3; ++j)
  {
    uint const * ev = topology(1, 0)(e[j]);
    dolfin_assert(ev);
    if (ev[0] != v[i] && ev[1] != v[i])
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge with index %d in triangle cell.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------
bool TriangleCell::check(Cell& cell) const
{
  bool ret = CellType::check(cell);

  // UFC convention: cell -> vertices in ascending order
  // These connectivities should always exist, catching assertion if it is not
  // the case is the right behaviour
  uint const * cell_verts = cell.entities(0);
  dolfin_assert(cell_verts);
  uint const num_cell_verts = this->num_vertices(this->dim());
  if(!is_sorted(cell_verts, cell_verts + num_cell_verts))
  {
    ret = false;
    warning("CellType::check : cell vertices are not in ascending order\n"
            "=> cell index = %d", cell.index());
  }

  // Check edge -> incident vertices mapping
  if (cell.mesh().topology().is_computed(1, 0))
  {
    uint const* v = cell.entities(0);
    dolfin_assert(v);
    uint const* e = cell.entities(1);
    dolfin_assert(e);
    MeshTopology const& topology = cell.mesh().topology();
    for (uint i = 0; i < 3; ++i)
    {
      uint const * ev = topology(1, 0)(e[i]);
      dolfin_assert(ev);
      for (uint j = 0; j < 2; ++j)
      {
        if (ev[j] != v[EIV[i][j]])
        {
          ret = false;
          warning("CellType : invalid edge -> incident vertices mapping");
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
