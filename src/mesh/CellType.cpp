// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/CellType.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/HexahedronCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexIterator.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
CellType::CellType(std::string const& name, CellType::Type cell_type,
                   CellType::Type facet_type) :
    name_(name),
    cell_type(cell_type),
    facet_type(facet_type),
    ufl_(CellType::ufldomain(cell_type))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CellType::~CellType()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Array<CellType*> CellType::create_all()
{
  Array<CellType *> ret;
  ret.push_back(CellType::create(CellType::interval));
  ret.push_back(CellType::create(CellType::triangle));
  ret.push_back(CellType::create(CellType::tetrahedron));
  ret.push_back(CellType::create(CellType::quadrilateral));
  ret.push_back(CellType::create(CellType::hexahedron));
  return ret;
}
//-----------------------------------------------------------------------------
Array<CellType*> CellType::create_simplex()
{
  Array<CellType *> ret;
  ret.push_back(CellType::create(CellType::interval));
  ret.push_back(CellType::create(CellType::triangle));
  ret.push_back(CellType::create(CellType::tetrahedron));
  return ret;
}
//-----------------------------------------------------------------------------
CellType* CellType::create_simplex(uint dim)
{
  switch (dim)
    {
    case 1:
      return CellType::create(CellType::interval);
      break;
    case 2:
      return CellType::create(CellType::triangle);
      break;
    case 3:
      return CellType::create(CellType::tetrahedron);
      break;
    default:
      error("Unknown simplex type for dimension: %d.", dim);
      break;
    }
  return NULL;
}
//-----------------------------------------------------------------------------
Array<CellType*> CellType::create_hypercube()
{
  Array<CellType *> ret;
  ret.push_back(CellType::create(CellType::interval));
  ret.push_back(CellType::create(CellType::quadrilateral));
  ret.push_back(CellType::create(CellType::hexahedron));
  return ret;
}
//-----------------------------------------------------------------------------
CellType* CellType::create_hypercube(uint dim)
{
  switch (dim)
    {
    case 1:
      return CellType::create(CellType::interval);
      break;
    case 2:
      return CellType::create(CellType::quadrilateral);
      break;
    case 3:
      return CellType::create(CellType::hexahedron);
      break;
    default:
      error("Unknown hypercube type for dimension: %d.", dim);
      break;
    }
  return NULL;
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

  return NULL;
}
//-----------------------------------------------------------------------------
CellType* CellType::create(std::string const& type)
{
  if (type == "point")
  {
    return new PointCell();
  }
  if (type == "interval")
  {
    return new IntervalCell();
  }
  else if (type == "triangle")
  {
    return new TriangleCell();
  }
  else if (type == "tetrahedron")
  {
    return new TetrahedronCell();
  }
  else if (type == "quadrilateral")
  {
    return new QuadrilateralCell();
  }
  else if (type == "hexahedron")
  {
    return new HexahedronCell();
  }
  else
  {
    error("Unknown cell type: \"%s\".", type.c_str());
  }

  return NULL;
}
//-----------------------------------------------------------------------------
CellType* CellType::create(ufl::Cell const& cell)
{
  switch (cell.domain().type())
    {
    case ufl::Domain::vertex:
      return new PointCell();
    case ufl::Domain::interval:
      return new IntervalCell();
    case ufl::Domain::triangle:
      return new TriangleCell();
    case ufl::Domain::tetrahedron:
      return new TetrahedronCell();
    case ufl::Domain::quadrilateral:
      return new QuadrilateralCell();
    case ufl::Domain::hexahedron:
      return new HexahedronCell();
    default:
      error("Unknown UFL domain type: %d.", type);
      break;
    }

  return NULL;
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
CellType::Type CellType::type(std::string const& type)
{
  if (type == "point")
  {
    return point;
  }
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

  return point;
}
//-----------------------------------------------------------------------------
std::string CellType::str(CellType::Type type)
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
std::string const& CellType::str() const
{
  return name_;
}
//-----------------------------------------------------------------------------
ufl::Domain::Type CellType::ufldomain(CellType::Type type)
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
bool CellType::check(Cell& cell) const
{
  // Throw a hard error
  if(cell.type() != this->cellType())
  {
    error("CellType::check : mismatch of cell type");
  }

  // UFC convention: edge -> vertices in ascending order
  if(cell.dim() < 2)
  {
    return true;
  }
  bool ret = true;
  if (cell.mesh().topology().connectivity(1, 0))
  {
    Array<uint> const & cell_edges = cell.entities(1);
    dolfin_assert( not cell_edges.empty() );
    uint const num_cell_edges = this->num_entities(1);
    for (uint e = 0; e < num_cell_edges; ++e)
    {
      Array<uint> const & edge_verts = cell.mesh().topology()(1, 0)[cell_edges[e]];
      dolfin_assert( not edge_verts.empty() );
      if (edge_verts[1] < edge_verts[0])
      {
        ret = false;
        warning("CellType::check : edge vertices are not in ascending order");
      }
    }
  }

  return ret;
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
