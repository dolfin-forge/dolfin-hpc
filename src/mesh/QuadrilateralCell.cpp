// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
//

#include <dolfin/mesh/QuadrilateralCell.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const QuadrilateralCell::NE[3][3] =
{ { 1, 0, 0 }, { 2, 1, 0 }, { 4, 4, 1 } };

// UFC: Vertex Coordinates
real const QuadrilateralCell::VC[4][2] =
{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

// UFC: Edge - Incident Vertices
uint const QuadrilateralCell::EIV[4][2] =
{ { 2, 3 }, { 1, 2 }, { 0, 3 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
uint const QuadrilateralCell::ENV[4][2] =
{ { 0, 1 }, { 0, 3 }, { 1, 2 }, { 2, 3 } };

//-----------------------------------------------------------------------------
QuadrilateralCell::QuadrilateralCell() :
    CellType("quadrilateral", CellType::quadrilateral, CellType::interval)
{
}
//-----------------------------------------------------------------------------
QuadrilateralCell::~QuadrilateralCell()
{
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::dim() const
{
  return 2;
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[2][dim];
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::num_entities(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD);
  dolfin_assert(d1 <= TD);
  return NE[d0][d1];
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim][0];
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get the coordinates of vertices v0, v1 and v2
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v2 = geometry.x(vertices[2]);

  // Check whether (v0v1, v0v2) is counter-clockwise
  return (
      ((v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0]))
          < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges
  if (dim != 1)
  {
    error("Invalid topological dimension for creation of entities: %d.", dim);
  }

  // Create the four edges using UFC ordering
  e[0][0] = v[2];
  e[0][1] = v[3];
  e[1][0] = v[1];
  e[1][1] = v[2];
  e[2][0] = v[0];
  e[2][1] = v[3];
  e[3][0] = v[0];
  e[3][1] = v[1];
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::order_entities(MeshTopology& topology, uint i) const
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
    for (uint i = 0; i < 4; ++i)
    {
      uint* edge_vertices = topology(1, 0)(cell_edges[i]);
      std::sort(edge_vertices, edge_vertices + 2);
    }
  }

  // Vertices are supposed to have be ordered appropriately to obtain a
  // conforming quadrilateral

  // Sort local edges on cell after non-incident vertex, connectivity 2 - 1
  if (topology.is_computed(2, 1))
  {
    dolfin_assert(topology.is_computed(2, 1));

    // Get cell vertices and edges
    uint* cell_vertices = topology(2, 0)(i);
    uint* cell_edges = topology(2, 1)(i);

    // Loop on non-incident vertices on cell as lexicographical pairs
    // (i, j): (0,1) (0,3) (1,2) (2,3)
    uint m = 0;
    for (uint t = 0; t < 4; ++t)
    {
      uint v0 = ENV[t][0];
      uint v1 = ENV[t][1];

      // Loop on edges
      for (uint k = m; k < 4; ++k)
      {
        // Get local vertices on edge
        uint* edge_vertices = topology(1, 0)(cell_edges[k]);

        // Check if the ith and jth vertex of the cell are non-incident on edge k
#if __SUNPRO_CC
        int n1 = 0;
        int n2 = 0;
        std::count(edge_vertices, edge_vertices+2, cell_vertices[v0], n1);
        std::count(edge_vertices, edge_vertices+2, cell_vertices[v1], n2);
        if (!n1 && !n2 )
#else
        if (!std::count(edge_vertices, edge_vertices + 2, cell_vertices[v0])
         && !std::count(edge_vertices, edge_vertices + 2, cell_vertices[v1]))
#endif
        {
          // Swap edge numbers
          uint tmp = cell_edges[m];
          cell_edges[m] = cell_edges[k];
          cell_edges[k] = tmp;
          ++m;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::order_facet(uint vertices[], Facet& facet) const
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
bool QuadrilateralCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  // Do not order cell - vertices connectivities
  return (d0 > 0 && d0 > d1) && !(d0 == TD && d1 == 0);
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::initialize_connectivities(Mesh& mesh) const
{
  mesh.init(1, 0);
  mesh.init(2, 0);
  mesh.init(2, 1);
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::refine_cell(Cell& cell, MeshEditor& editor,
                                   uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  uint const * v = cell.entities(0);
  dolfin_assert(v);
  uint const * e = cell.entities(1);
  dolfin_assert(e);

  // Compute indices for the nine new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const v3 = v[3];
  uint const eoffset = cell.mesh().size(0);
  uint const e0 = eoffset + e[findEdge(0, cell)];
  uint const e1 = eoffset + e[findEdge(1, cell)];
  uint const e2 = eoffset + e[findEdge(2, cell)];
  uint const e3 = eoffset + e[findEdge(3, cell)];
  uint const coffset = eoffset + cell.mesh().size(1);
  uint const c0 = coffset + cell.index();

  // Add the four new cells
  uint cv0[4] = { v0, e3, c0, e2 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint cv1[4] = { e3, v1, e1, c0};
  editor.add_cell(current_cell++, &cv1[0]);
  uint cv2[4] = { c0, e1, v2, e0 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint cv3[4] = { e2, c0, e0, v3 };
  editor.add_cell(current_cell++, &cv3[0]);
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::num_refined_cells() const
{
  return 4;
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return 1;
}
//-----------------------------------------------------------------------------
real QuadrilateralCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[2][0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * a = geometry.x(vertices[0]);
  real const * b = geometry.x(vertices[1]);
  real const * c = geometry.x(vertices[2]);
  real const * d = geometry.x(vertices[3]);

  switch (geometry.dim())
    {
    case 2:
    {
      return 0.5 *
          std::fabs((c[0] - a[0]) * (d[1] - b[1])
                    - (c[1] - a[1]) * (d[0] - b[0]));
    }
    break;
    case 3:
      {
        real c0 = (c[1] - a[1]) * (d[2] - b[2]) - (c[2] - a[2]) * (d[1] - b[1]);
        real c1 = (c[2] - a[2]) * (d[0] - b[0]) - (c[0] - a[0]) * (d[2] - b[2]);
        real c2 = (c[0] - a[0]) * (d[1] - b[1]) - (c[1] - a[1]) * (d[0] - b[0]);
        return 0.5 * std::sqrt(c0 * c0 + c1 * c1 + c2 * c2);
      }
      break;
    }

  return 0.0;
}
//-----------------------------------------------------------------------------
real QuadrilateralCell::diameter(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[2][0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);

  // Compute maximum diagonal length
  real d0 = 0.0;
  real d1 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    d0 += (x3[i] - x0[i]) * (x3[i] - x0[i]);
    d1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
  }
  return std::sqrt(std::max(d0, d1));
}
//-----------------------------------------------------------------------------
real QuadrilateralCell::circumradius(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[2][0]);

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);

  // Compute circumradius assuming that the quadrilateral is cyclic
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real d = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    a += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    b += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    c += (x3[i] - x2[i]) * (x3[i] - x2[i]);
    d += (x0[i] - x3[i]) * (x0[i] - x3[i]);
  }
  a = std::sqrt(a);
  b = std::sqrt(b);
  c = std::sqrt(c);
  d = std::sqrt(d);

  // Using Brahmagupta's formula for the volume
  real const s = 0.5 * (a + b + c + d);
  return 0.25 *
      std::sqrt((a*c+b*d)*(a*d+b*c)*(a*b+c*d)/((s-a)*(s-b)*(s-c)*(s-d)));
}
//-----------------------------------------------------------------------------
real QuadrilateralCell::inradius(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[2][0]);

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);

  // Compute inradius assuming that the quadrilateral is tangent
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real d = 0.0;
  real p = 0.0;
  real q = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    a += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    b += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    c += (x3[i] - x2[i]) * (x3[i] - x2[i]);
    d += (x0[i] - x3[i]) * (x0[i] - x3[i]);
    p += (x2[i] - x0[i]) * (x2[i] - x0[i]);
    q += (x3[i] - x1[i]) * (x3[i] - x1[i]);
  }
  a = std::sqrt(a);
  b = std::sqrt(b);
  c = std::sqrt(c);
  d = std::sqrt(d);
  p = std::sqrt(p);
  q = std::sqrt(q);

  return 0.5
      * std::sqrt(
          4.0 * p * p * q * q
              - (a * a - b * b + c * c - d * d)
                  * (a * a - b * b + c * c - d * d)) / (a + b + c + d);
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::midpoint(MeshEntity const& entity, real * p) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[2][0]);

  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  real const* x3 = geometry.x(vertices[3]);
  uint const gdim = geometry.dim();
  for (uint d = 0; d < gdim; ++d)
  {
    p[d] = 0.25 * (x0[d] + x1[d] + x2[d] + x3[d]);
  }
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::normal(Cell const& cell, uint facet, real * n) const
{
  dolfin_assert(cell.type() == this->cell_type);

  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();
  // Get coordinates of first non-incident vertex
  real const * p0 = geometry.x(c.entities(0)[ENV[facet][0]]);
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
      error("QuadrilateralCell : facet normal not implemented in R^%u", gdim);
      break;
  }
}
//-----------------------------------------------------------------------------
real QuadrilateralCell::facet_area(Cell const& cell, uint facet) const
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
bool QuadrilateralCell::intersects(MeshEntity const& e, Point const& p) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[2][0]);

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
   */

  error("Collision of quadrilateral with point no implemented.");

  return true;
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::intersects(MeshEntity const& e, Point const& p1,
                                   Point const& p2) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[2][0]);

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
   */

  error("Collision of quadrilateral with segment not implemented");

  return true;
}
//-----------------------------------------------------------------------------
std::string QuadrilateralCell::description() const
{
  return std::string("quadrilateral (hypercube of topological dimension 2)");
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::quadrilateral, 2, DOLFIN_COMM_SELF);
  me.init_vertices(4);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.add_vertex(2, VC[2]);
  me.add_vertex(3, VC[3]);
  me.init_cells(1);
  uint const cv0[4] = { 0, 1, 2, 3 };
  me.add_cell(0, cv0);
  me.close();
}
//-----------------------------------------------------------------------------
real const * QuadrilateralCell::reference_vertex(uint i) const
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::disp() const
{
  message("QuadrilateralCell");
  begin(  "-----------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::check(Cell& cell) const
{
  bool ret = CellType::check(cell);

  // Check edge -> incident vertices mapping
  if (cell.mesh().topology().is_computed(1, 0))
  {
    uint const* v = cell.entities(0);
    dolfin_assert(v);
    uint const* e = cell.entities(1);
    dolfin_assert(e);
    MeshTopology const& topology = cell.mesh().topology();
    for (uint i = 0; i < 4; ++i)
    {
      uint const * ev = topology(1, 0)(e[i]);
      dolfin_assert(ev);
      for (uint j = 0; j < 2; ++j)
      {
        if (ev[j] != v[EIV[i][0]] && ev[j] != v[EIV[i][1]])
        {
          ret = false;
          warning("CellType : invalid edge -> incident vertices mapping\n"
                "e[%u] = (v%u, v%u)", i, ev[0], ev[1]);
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::findEdge(uint i, Cell const& cell) const
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  uint const * v = cell.entities(0);
  dolfin_assert(v);
  uint const * e = cell.entities(1);
  dolfin_assert(e);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  uint const v0 = v[ENV[i][0]];
  uint const v1 = v[ENV[i][1]];
  for (uint j = 0; j < 4; ++j)
  {
    uint const * ev = topology(1, 0)(e[j]);
    dolfin_assert(ev);
    if (ev[0] != v0 && ev[0] != v1 && ev[1] != v0 && ev[1] != v1)
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge with index %d in quadrilateral.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
