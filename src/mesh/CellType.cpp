// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-05
// Last changed: 2006-10-16

#include <dolfin/mesh/CellType.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>

#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/HexahedronCell.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
CellType::CellType(CellType::Type cell_type, CellType::Type facet_type) :
    cell_type(cell_type),
    facet_type(facet_type),
    ufl_(CellType::type2ufldomain(cell_type))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CellType::~CellType()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CellType* CellType::create(CellType::Type type)
{
  switch (type)
    {
    case point:
      return new PointCell();
    case interval:
      return new IntervalCell();
    case triangle:
      return new TriangleCell();
    case tetrahedron:
      return new TetrahedronCell();
    case quadrilateral:
      return new QuadrilateralCell();
    case hexahedron:
      return new HexahedronCell();
    default:
      error("Unknown cell type: %d.", type);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
CellType* CellType::create(std::string type)
{
  return create(string2type(type));
}
//-----------------------------------------------------------------------------
CellType::Type CellType::string2type(std::string type)
{
  if (type == "interval")
  {
    return interval;
  }
  else if (type == "triangle")
  {
    return triangle;
  }
  else if (type == "tetrahedron")
  {
    return tetrahedron;
  }
  else if (type == "quadrilateral")
  {
      return quadrilateral;
  }
  else if (type == "hexahedron")
  {
    return hexahedron;
  }
  else
  {
    error("Unknown cell type: \"%s\".", type.c_str());
  }

  return interval;
}
//-----------------------------------------------------------------------------
bool CellType::intersects(MeshEntity& entity, Cell& c) const
{
  for (VertexIterator vi(entity); !vi.end(); ++vi)
  {
    Point p = vi->point();

    if (intersects(c, p))
    {
      return true;
    }
  }

  for (VertexIterator vi(c); !vi.end(); ++vi)
  {
    Point p = vi->point();

    if (intersects(entity, p))
    {
      return true;
    }
  }

  return false;
}
//-----------------------------------------------------------------------------
std::string CellType::type2string(CellType::Type type)
{
  switch (type)
    {
    case point:
      return "point";
    case interval:
      return "interval";
    case triangle:
      return "triangle";
    case tetrahedron:
      return "tetrahedron";
    case quadrilateral:
      return "quadrilateral";
    case hexahedron:
      return "hexahedron";
    default:
      error("Unknown cell type: %d.", type);
      break;
    }

  return "";
}
//-----------------------------------------------------------------------------
ufl::Domain::Type CellType::type2ufldomain(CellType::Type type)
{
  switch (type)
    {
    case CellType::point:
      return ufl::Domain::vertex;
    case CellType::interval:
      return ufl::Domain::interval;
    case CellType::triangle:
      return ufl::Domain::triangle;
    case CellType::tetrahedron:
      return ufl::Domain::tetrahedron;
    case CellType::quadrilateral:
      return ufl::Domain::quadrilateral;
    case CellType::hexahedron:
      return ufl::Domain::hexahedron;
    default:
      error("Unknown cell type: %d.", type);
      break;
    }

  return ufl::Domain::None;
}
//-----------------------------------------------------------------------------
void CellType::check(Cell& cell) const
{
  if(cell.type() != this->cellType())
  {
    error("CellType::check : mismatch of cell type");
  }

  // UFC convention: cell -> vertices in ascending order
  uint const * cell_verts = cell.entities(0);
  dolfin_assert(cell_verts);
  uint const num_cell_verts = this->numVertices(this->dim());
  if(!is_sorted(cell_verts, cell_verts + num_cell_verts))
  {
    error("CellType::check : cell vertices are not in ascending order\n"
          "=> cell index = %d", cell.index());
  }

  // UFC convention: edge -> vertices in ascending order
  if(cell.dim() < 2)
  {
    return;
  }
  uint const * cell_edges = cell.entities(1);
  dolfin_assert(cell_edges);
  uint const num_cell_edges = this->numEntities(1);
  uint const num_edge_verts = this->numVertices(1);
  for (uint e = 0; e < num_cell_edges; ++e)
  {
    uint const * edge_verts = cell.mesh().topology()(1, 0)(cell_edges[e]);
    dolfin_assert(edge_verts);
    if (edge_verts[1] < edge_verts[0])
    {
      error("CellType::check : edge vertices are not in ascending order");
    }
  }
}
//-----------------------------------------------------------------------------
uint const * CellType::is_sorted_until(uint const * begin, uint const * end)
{
  if (begin == end)
  {
    return begin;
  }
  uint const * next = begin;
  while (++next != end)
  {
    if (*next < *begin)
    {
      return next;
    }
    ++begin;
  }
  return end;
}
//-----------------------------------------------------------------------------
bool CellType::is_sorted(uint const * begin, uint const * end)
{
  return (is_sorted_until(begin, end) == end);
}
//-----------------------------------------------------------------------------
bool CellType::pattern_applies(Cell& cell)const
{
  return (cell.type() == this->cellType());
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
