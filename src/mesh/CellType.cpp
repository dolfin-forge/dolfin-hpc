// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/CellType.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/HexahedronCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

#include <ufc.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------

CellType::CellType( std::string const & name,
                    CellType::Type      cell_type,
                    CellType::Type      facet_type )
  : name_( name )
  , cell_type( cell_type )
  , facet_type( facet_type )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

CellType::~CellType()
{
  // Do nothing
}

//-----------------------------------------------------------------------------

auto CellType::create_all() -> std::vector< CellType * >
{
  std::vector< CellType * > ret;
  ret.push_back( CellType::create( CellType::interval ) );
  ret.push_back( CellType::create( CellType::triangle ) );
  ret.push_back( CellType::create( CellType::tetrahedron ) );
  ret.push_back( CellType::create( CellType::quadrilateral ) );
  ret.push_back( CellType::create( CellType::hexahedron ) );
  return ret;
}

//-----------------------------------------------------------------------------

auto CellType::create_simplex() -> std::vector< CellType * >
{
  std::vector< CellType * > ret;
  ret.push_back( CellType::create( CellType::interval ) );
  ret.push_back( CellType::create( CellType::triangle ) );
  ret.push_back( CellType::create( CellType::tetrahedron ) );
  return ret;
}

//-----------------------------------------------------------------------------

auto CellType::create_simplex( size_t dim ) -> CellType *
{
  switch ( dim )
  {
    case 1:
      return CellType::create( CellType::interval );
      break;
    case 2:
      return CellType::create( CellType::triangle );
      break;
    case 3:
      return CellType::create( CellType::tetrahedron );
      break;
    default:
      error( "Unknown simplex type for dimension: %d.", dim );
      break;
  }
  return nullptr;
}

//-----------------------------------------------------------------------------

auto CellType::create_hypercube() -> std::vector< CellType * >
{
  std::vector< CellType * > ret;
  ret.push_back( CellType::create( CellType::interval ) );
  ret.push_back( CellType::create( CellType::quadrilateral ) );
  ret.push_back( CellType::create( CellType::hexahedron ) );
  return ret;
}

//-----------------------------------------------------------------------------

auto CellType::create_hypercube( size_t dim ) -> CellType *
{
  switch ( dim )
  {
    case 1:
      return CellType::create( CellType::interval );
      break;
    case 2:
      return CellType::create( CellType::quadrilateral );
      break;
    case 3:
      return CellType::create( CellType::hexahedron );
      break;
    default:
      error( "Unknown hypercube type for dimension: %d.", dim );
      break;
  }
  return nullptr;
}

//-----------------------------------------------------------------------------

auto CellType::create( CellType::Type type ) -> CellType *
{
  switch ( type )
  {
    case point:
      return new PointCell();
      break;
    case interval:
      return new IntervalCell();
      break;
    case triangle:
      return new TriangleCell();
      break;
    case tetrahedron:
      return new TetrahedronCell();
      break;
    case quadrilateral:
      return new QuadrilateralCell();
      break;
    case hexahedron:
      return new HexahedronCell();
      break;
    default:
      error( "Unknown cell type: %d.", type );
      break;
  }

  return nullptr;
}

//-----------------------------------------------------------------------------

auto CellType::create( ufc::shape type ) -> CellType *
{
  switch ( type )
  {
    case ufc::shape::vertex:
      return new PointCell();
      break;
    case ufc::shape::interval:
      return new IntervalCell();
      break;
    case ufc::shape::triangle:
      return new TriangleCell();
      break;
    case ufc::shape::tetrahedron:
      return new TetrahedronCell();
      break;
    case ufc::shape::quadrilateral:
      return new QuadrilateralCell();
      break;
    case ufc::shape::hexahedron:
      return new HexahedronCell();
      break;
    default:
      error( "Unknown cell type: %d.", type );
      break;
  }

  return nullptr;
}

//-----------------------------------------------------------------------------

auto CellType::create( std::string const & type ) -> CellType *
{
  if ( type == "point" )
  {
    return new PointCell();
  }
  if ( type == "interval" )
  {
    return new IntervalCell();
  }
  else if ( type == "triangle" )
  {
    return new TriangleCell();
  }
  else if ( type == "tetrahedron" )
  {
    return new TetrahedronCell();
  }
  else if ( type == "quadrilateral" )
  {
    return new QuadrilateralCell();
  }
  else if ( type == "hexahedron" )
  {
    return new HexahedronCell();
  }
  else
  {
    error( "Unknown cell type: \"%s\".", type.c_str() );
  }

  return nullptr;
}

//-----------------------------------------------------------------------------

auto CellType::intersects( MeshEntity & entity, Cell & c ) const -> bool
{
  for ( VertexIterator vi( entity ); !vi.end(); ++vi )
  {
    Point p = vi->point();

    if ( intersects( c, p ) )
    {
      return true;
    }
  }

  for ( VertexIterator vi( c ); !vi.end(); ++vi )
  {
    Point p = vi->point();

    if ( intersects( entity, p ) )
    {
      return true;
    }
  }

  return false;
}

//-----------------------------------------------------------------------------

auto CellType::type( std::string const & type ) -> CellType::Type
{
  if ( type == "point" )
  {
    return point;
  }
  if ( type == "interval" )
  {
    return interval;
  }
  else if ( type == "triangle" )
  {
    return triangle;
  }
  else if ( type == "tetrahedron" )
  {
    return tetrahedron;
  }
  else if ( type == "quadrilateral" )
  {
    return quadrilateral;
  }
  else if ( type == "hexahedron" )
  {
    return hexahedron;
  }
  else
  {
    error( "Unknown cell type: \"%s\".", type.c_str() );
  }

  return point;
}

//-----------------------------------------------------------------------------

auto CellType::str( CellType::Type type ) -> std::string
{
  switch ( type )
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
      error( "Unknown cell type: %d.", type );
      break;
  }

  return "";
}

//-----------------------------------------------------------------------------

auto CellType::str() const -> std::string const &
{
  return name_;
}

//-----------------------------------------------------------------------------

auto CellType::check( Cell & cell ) const -> bool
{
  // Throw a hard error
  if ( cell.type() != this->cellType() )
  {
    error( "CellType::check : mismatch of cell type" );
  }

  // UFC convention: edge -> vertices in ascending order
  if ( cell.dim() < 2 )
  {
    return true;
  }
  bool ret = true;
  if ( cell.mesh().topology().connectivity( 1, 0 ) )
  {
    std::vector< size_t > const & cell_edges = cell.entities( 1 );
    dolfin_assert( not cell_edges.empty() );
    size_t const num_cell_edges = this->num_entities( 1 );
    for ( size_t e = 0; e < num_cell_edges; ++e )
    {
      std::vector< size_t > const & edge_verts =
        cell.mesh().topology()( 1, 0 )[cell_edges[e]];
      dolfin_assert( not edge_verts.empty() );
      if ( edge_verts[1] < edge_verts[0] )
      {
        ret = false;
        warning( "CellType::check : edge vertices are not in ascending order" );
      }
    }
  }

  return ret;
}

//-----------------------------------------------------------------------------

auto CellType::is_sorted_until( size_t const * begin, size_t const * end )
  -> size_t const *
{
  if ( begin == end )
  {
    return begin;
  }
  size_t const * next = begin;
  while ( ++next != end )
  {
    if ( *next < *begin )
    {
      return next;
    }
    ++begin;
  }
  return end;
}

//-----------------------------------------------------------------------------

auto CellType::is_sorted( size_t const * begin, size_t const * end ) -> bool
{
  return ( is_sorted_until( begin, end ) == end );
}

//-----------------------------------------------------------------------------

auto CellType::pattern_applies( Cell & cell ) const -> bool
{
  return ( cell.type() == this->cellType() );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
