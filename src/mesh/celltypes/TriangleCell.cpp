// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/celltypes/TriangleCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/GeometricPredicates.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
size_t const TriangleCell::NE[3][3] = { { 1, 0, 0 }, { 2, 1, 0 }, { 3, 3, 1 } };

// UFC: Vertex Coordinates
real const TriangleCell::VC[3][2] = { { 0.0, 0.0 },
                                      { 1.0, 0.0 },
                                      { 0.0, 1.0 } };

// UFC: Edge - Incident Vertices
size_t const TriangleCell::EIV[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
size_t const TriangleCell::ENV[3][1] = { { 0 }, { 1 }, { 2 } };

//-----------------------------------------------------------------------------
TriangleCell::TriangleCell()
  : CellType( "triangle", CellType::triangle, CellType::interval )
{
}
//-----------------------------------------------------------------------------
void TriangleCell::create_entities( size_t **      e,
                                    size_t         dim,
                                    size_t const * v ) const
{
  // We only need to know how to create edges
  if ( dim != 1 )
  {
    error( "Invalid topological dimension for creation of entities: %d.", dim );
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
void TriangleCell::order_entities( MeshTopology & topology, size_t i ) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1
  dolfin_assert( topology.type().cellType() == this->cell_type );

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if ( topology.connectivity( 1, 0 ) )
  {
    dolfin_assert( topology.connectivity( 2, 1 ) );

    // Get edges
    Array< size_t > const & cell_edges = topology( 2, 1 )[i];

    // Sort vertices on each edge
    for ( size_t i = 0; i < 3; ++i )
    {
      Array< size_t > & edge_vertices = topology( 1, 0 )[cell_edges[i]];
      std::sort( edge_vertices.data(), edge_vertices.data() + 2 );
    }
  }

  // Sort local vertices on cell in ascending order, connectivity 2 - 0
  if ( topology.connectivity( 2, 0 ) )
  {
    Array< size_t > & cell_vertices = topology( 2, 0 )[i];
    std::sort( cell_vertices.data(), cell_vertices.data() + 3 );
  }

  // Sort local edges on cell after non-incident vertex, connectivity 2 - 1
  if ( topology.connectivity( 2, 1 ) )
  {
    dolfin_assert( topology.connectivity( 2, 1 ) );

    // Get cell vertices and edges
    Array< size_t > const & cell_vertices = topology( 2, 0 )[i];
    Array< size_t > &       cell_edges    = topology( 2, 1 )[i];

    // Loop over vertices on cell
    for ( size_t i = 0; i < 3; ++i )
    {
      // Loop over edges on cell
      for ( size_t j = i; j < 3; ++j )
      {
        Array< size_t > const & edge_vertices = topology( 1, 0 )[cell_edges[j]];

        // Check if the ith vertex of the cell is non-incident with edge j
        if ( std::count( edge_vertices.data(),
                         edge_vertices.data() + 2,
                         cell_vertices[i] )
             == 0 )
        {
          // Swap edge numbers
          size_t tmp    = cell_edges[i];
          cell_edges[i] = cell_edges[j];
          cell_edges[j] = tmp;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void TriangleCell::order_facet( size_t vertices[], Facet & facet ) const
{
  // Get mesh
  Mesh & mesh = facet.mesh();

  // Get the vertex opposite to the facet (the one we remove)
  size_t     vertex = 0;
  const Cell cell( mesh, facet.entities( mesh.topology_dimension() )[0] );
  for ( size_t i = 0; i < cell.num_entities( 0 ); i++ )
  {
    bool not_in_facet = true;
    vertex            = cell.entities( 0 )[i];
    for ( size_t j = 0; j < facet.num_entities( 0 ); j++ )
    {
      if ( vertex == facet.entities( 0 )[j] )
      {
        not_in_facet = false;
        break;
      }
    }
    if ( not_in_facet )
    {
      break;
    }
  }

  // Order
  Point const p  = mesh.geometry().point( vertex );
  Point       p0 = mesh.geometry().point( facet.entities( 0 )[0] );
  Point       p1 = mesh.geometry().point( facet.entities( 0 )[1] );
  Point       v  = p1 - p0;
  Point       n( v[1], -v[0], 0.0 );
  if ( n.dot( p0 - p ) < 0.0 )
  {
    std::swap( vertices[0], vertices[1] );
  }
}
//-----------------------------------------------------------------------------
void TriangleCell::refine_cell( Cell &       cell,
                                MeshEditor & editor,
                                size_t &     current_cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get vertices and edges
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( !v.empty() );
  Array< size_t > const & e = cell.entities( 1 );
  dolfin_assert( !e.empty() );

  // Compute indices for the six new vertices
  size_t const v0     = v[0];
  size_t const v1     = v[1];
  size_t const v2     = v[2];
  size_t const offset = cell.mesh().size( 0 );
  size_t const e0     = offset + e[findEdge( 0, cell )];
  size_t const e1     = offset + e[findEdge( 1, cell )];
  size_t const e2     = offset + e[findEdge( 2, cell )];

  // Add the four new cells
  size_t cv0[3] = { v0, e2, e1 };
  editor.add_cell( current_cell++, &cv0[0] );
  size_t cv1[3] = { v1, e0, e2 };
  editor.add_cell( current_cell++, &cv1[0] );
  size_t cv2[3] = { v2, e1, e0 };
  editor.add_cell( current_cell++, &cv2[0] );
  size_t cv3[3] = { e0, e1, e2 };
  editor.add_cell( current_cell++, &cv3[0] );
}
//-----------------------------------------------------------------------------
void TriangleCell::normal( Cell const & cell, size_t facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  Cell &               c = const_cast< Cell & >( cell );
  Facet                f( c.mesh(), c.entities( 1 )[facet] );
  MeshGeometry const & geometry = cell.mesh().geometry();
  // Get coordinates of opposite vertex
  real const * p0 = geometry.x( c.entities( 0 )[facet] );
  // Get coordinates of edge vertices
  Array< size_t > const & vertices = f.entities( 0 );
  real const *            p1       = geometry.x( vertices[0] );
  real const *            p2       = geometry.x( vertices[1] );
  size_t const            gdim     = geometry.dim();
  switch ( gdim )
  {
    case 2:
    {
      n[0]          = p2[1] - p1[1];
      n[1]          = p1[0] - p2[0];
      real const nn = std::sqrt( n[0] * n[0] + n[1] * n[1] );
      n[0] /= nn;
      n[1] /= nn;
      // Flip direction of normal so it points outward
      if ( ( n[0] * ( p1[0] - p0[0] ) + n[1] * ( p1[1] - p0[1] ) ) < 0.0 )
      {
        n[0] *= -1.0;
        n[1] *= -1.0;
      }
    }
    break;
    case 3:
    {
      real ac = ( p2[0] - p1[0] ) * ( p2[0] - p0[0] )
                + ( p2[1] - p1[1] ) * ( p2[1] - p0[1] )
                + ( p2[2] - p1[2] ) * ( p2[2] - p0[2] );
      real ab = ( p2[0] - p1[0] ) * ( p1[0] - p0[0] )
                + ( p2[1] - p1[1] ) * ( p1[1] - p0[1] )
                + ( p2[2] - p1[2] ) * ( p1[2] - p0[2] );
      n[0]    = ( p1[0] - p0[0] ) * ac - ( p2[0] - p0[0] ) * ab;
      n[1]    = ( p1[1] - p0[1] ) * ac - ( p2[1] - p0[1] ) * ab;
      n[2]    = ( p1[2] - p0[2] ) * ac - ( p2[2] - p0[2] ) * ab;
      real nn = std::sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
      n[0] /= nn;
      n[1] /= nn;
      n[2] /= nn;
      dolfin_assert( ( n[0] * ( p1[0] - p0[0] ) + n[1] * ( p1[1] - p0[1] )
                       + n[2] * ( p1[2] - p0[2] ) )
                     > 0.0 );
    }
    break;
    default:
      error( "TriangleCell : facet normal not implemented in R^%u", gdim );
      break;
  }
}
//-----------------------------------------------------------------------------
auto TriangleCell::facet_area( Cell const & cell, size_t facet ) const -> real
{
  dolfin_assert( cell.type() == this->cell_type );

  Cell &                  c = const_cast< Cell & >( cell );
  Facet                   f( c.mesh(), c.entities( 1 )[facet] );
  MeshGeometry const &    geometry = cell.mesh().geometry();
  Array< size_t > const & vertices = f.entities( 0 );
  real const *            p0       = geometry.x( vertices[0] );
  real const *            p1       = geometry.x( vertices[1] );
  // Compute distance between vertices
  real         meas = 0.0;
  size_t const gdim = geometry.dim();
  for ( size_t d = 0; d < gdim; ++d )
  {
    meas += ( p1[d] - p0[d] ) * ( p1[d] - p0[d] );
  }
  return std::sqrt( meas );
}
//-----------------------------------------------------------------------------
auto TriangleCell::intersects( MeshEntity const & e, Point const & p ) const
  -> bool
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[2][0] );

  // Get mesh geometry
  MeshGeometry const & geometry = e.mesh().geometry();

  // Get global index of vertices of the triangle
  size_t const ort = orientation( ( Cell & ) e );
  size_t const v0  = e.entities( 0 )[0];
  size_t const v1  = ( ort == 0 ? e.entities( 0 )[1] : e.entities( 0 )[2] );
  size_t const v2  = ( ort == 0 ? e.entities( 0 )[2] : e.entities( 0 )[1] );

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x( v0 );
  real const * x1 = geometry.x( v1 );
  real const * x2 = geometry.x( v2 );

  // The return value of orient2d is defined as twice the signed measure of
  // the triangle defined by the points given as arguments.
  real tol = geometry.abs_tolerance( 2 );

  if ( geometry.dim() == 2 )
  {

    real l0 = std::sqrt( ( x1[0] - x0[0] ) * ( x1[0] - x0[0] )
                         + ( x1[1] - x0[1] ) * ( x1[1] - x0[1] ) );
    real l1 = std::sqrt( ( x2[0] - x1[0] ) * ( x2[0] - x1[0] )
                         + ( x2[1] - x1[1] ) * ( x2[1] - x1[1] ) );
    real l2 = std::sqrt( ( x0[0] - x2[0] ) * ( x0[0] - x2[0] )
                         + ( x0[1] - x2[1] ) * ( x0[1] - x2[1] ) );

    // Test orientation of p w.r.t. each edge
    real d1 = orient2d( x0, x1, &p[0] );
    if ( d1 < ( -2.0 * tol * l0 ) )
      return false;
    real d2 = orient2d( x1, x2, &p[0] );
    if ( d2 < ( -2.0 * tol * l1 ) )
      return false;
    real d3 = orient2d( x2, x0, &p[0] );
    if ( d3 < ( -2.0 * tol * l2 ) )
      return false;
  }
  else if ( geometry.dim() == 3 )
  {
    // real l0 = std::sqrt((x1[0] - x0[0])*(x1[0] - x0[0])
    //                     + (x1[1] - x0[1])*(x1[1] - x0[1])
    //                     + (x1[2] - x0[2])*(x1[2] - x0[2]));
    // real l1 = std::sqrt((x2[0] - x1[0])*(x2[0] - x1[0])
    //                     + (x2[1] - x1[1])*(x2[1] - x1[1])
    //                     + (x2[2] - x1[2])*(x2[2] - x1[2]));
    // real l2 = std::sqrt((x0[0] - x2[0])*(x0[0] - x2[0])
    //                     + (x0[1] - x2[1])*(x0[1] - x2[1])
    //                     + (x0[2] - x2[2])*(x0[2] - x2[2]));

    // real n0 = (x1[1] - x0[1])*(x2[2] - x0[2]) - (x1[2] - x0[2])*(x2[1] -
    // x0[1]); real n1 = (x1[2] - x0[2])*(x2[0] - x0[0]) - (x1[0] -
    // x0[0])*(x2[2] - x0[2]); real n2 = (x1[0] - x0[0])*(x2[1] - x0[1]) -
    // (x1[1] - x0[1])*(x2[0] - x0[0]); real zz = (n2*n2) / (n0*n0 + n1*n1 +
    // n2*n2);

    // Check direction of the normal to the triangle; |n.ez|/||n|| >= sqrt(2)/2
    real d1 = orient3d( x0, x1, x2, &p[0] );
    if ( std::abs( d1 ) > tol )
      return false;
  }
  else
  {
    error( "Collision of triangle with point only implemented in R^2 and R^3" );
  }

  return true;
}
//-----------------------------------------------------------------------------
auto TriangleCell::intersects( MeshEntity const & e,
                               Point const &      p1,
                               Point const &      p2 ) const -> bool
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[2][0] );

  // Get mesh geometry
  MeshGeometry const & geometry = e.mesh().geometry();

  // Get global index of vertices of the triangle
  size_t const ort = orientation( ( Cell & ) e );
  size_t const v0  = e.entities( 0 )[0];
  size_t const v1  = ( ort == 1 ? e.entities( 0 )[2] : e.entities( 0 )[1] );
  size_t const v2  = ( ort == 1 ? e.entities( 0 )[1] : e.entities( 0 )[2] );

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x( v0 );
  real const * x1 = geometry.x( v1 );
  real const * x2 = geometry.x( v2 );

  real d1, d2, d3;

  if ( geometry.dim() == 2 )
  {
    // Test orientation of each vertex w.r.t. p1-p2
    d1 = orient2d( &p1[0], &p2[0], x0 );
    d2 = orient2d( &p1[0], &p2[0], x1 );
    d3 = orient2d( &p1[0], &p2[0], x2 );

    if ( d1 < 0 && d2 < 0 && d3 < 0 )
      return false;
    if ( d1 > 0 && d2 > 0 && d3 > 0 )
      return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x0-x1:
    d1 = orient2d( x0, x1, &p1[0] );
    d2 = orient2d( x0, x1, &p2[0] );

    if ( d1 < 0 && d2 < 0 )
      return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x1-x2:
    d1 = orient2d( x1, x2, &p1[0] );
    d2 = orient2d( x1, x2, &p2[0] );

    if ( d1 < 0 && d2 < 0 )
      return false;

    // Line p1-p2 intersects triangle but both p1 and p2 are
    // on the negative side of x2-x0:
    d1 = orient2d( x2, x0, &p1[0] );
    d2 = orient2d( x2, x0, &p2[0] );

    if ( d1 < 0 && d2 < 0 )
      return false;
  }
  else
  {
    error( "Collision of triangle with segment only implemented in R^2" );
  }

  return true;
}
//-----------------------------------------------------------------------------
void TriangleCell::create_reference_cell( Mesh & mesh ) const
{
  MeshEditor me( mesh, CellType::triangle, 2, DOLFIN_COMM_SELF );
  me.init_vertices( 3 );
  me.add_vertex( 0, VC[0] );
  me.add_vertex( 1, VC[1] );
  me.add_vertex( 2, VC[2] );
  me.init_cells( 1 );
  size_t const cv0[3] = { 0, 1, 2 };
  me.add_cell( 0, cv0 );
  me.close();
}
//-----------------------------------------------------------------------------
auto TriangleCell::reference_vertex( size_t i ) const -> real const *
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------
auto TriangleCell::description() const -> std::string
{
  return std::string( "triangle (simplex of topological dimension 2)" );
}
//-----------------------------------------------------------------------------
void TriangleCell::disp() const
{
  message( "TriangleCell" );
  begin( "------------" );
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
auto TriangleCell::findEdge( size_t i, Cell const & cell ) const -> size_t
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( not v.empty() );
  Array< size_t > const & e = cell.entities( 1 );
  dolfin_assert( not e.empty() );

  // Look for edge satisfying ordering convention
  MeshTopology const & topology = cell.mesh().topology();
  for ( size_t j = 0; j < 3; ++j )
  {
    Array< size_t > const & ev = topology( 1, 0 )[e[j]];
    dolfin_assert( !ev.empty() );
    if ( ev[0] != v[i] && ev[1] != v[i] )
    {
      return j;
    }
  }

  // We should not reach this
  error( "Unable to find edge with index %d in triangle cell.", cell.index() );

  return 0;
}
//-----------------------------------------------------------------------------
auto TriangleCell::check( Cell & cell ) const -> bool
{
  bool ret = CellType::check( cell );

  // UFC convention: cell -> vertices in ascending order
  // These connectivities should always exist, catching assertion if it is not
  // the case is the right behaviour
  Array< size_t > const & cell_verts = cell.entities( 0 );
  dolfin_assert( not cell_verts.empty() );
  size_t const num_cell_verts = this->num_vertices( this->dim() );
  if ( !is_sorted( cell_verts.data(), cell_verts.data() + num_cell_verts ) )
  {
    ret = false;
    warning( "CellType::check : cell vertices are not in ascending order\n"
             "=> cell index = %d",
             cell.index() );
  }

  // Check edge -> incident vertices mapping
  if ( cell.mesh().topology().connectivity( 1, 0 ) )
  {
    Array< size_t > const & v = cell.entities( 0 );
    dolfin_assert( !v.empty() );
    Array< size_t > const & e = cell.entities( 1 );
    dolfin_assert( !e.empty() );
    MeshTopology const & topology = cell.mesh().topology();
    for ( size_t i = 0; i < 3; ++i )
    {
      Array< size_t > const & ev = topology( 1, 0 )[e[i]];
      dolfin_assert( not ev.empty() );
      for ( size_t j = 0; j < 2; ++j )
      {
        if ( ev[j] != v[EIV[i][j]] )
        {
          ret = false;
          warning( "CellType : invalid edge -> incident vertices mapping" );
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
