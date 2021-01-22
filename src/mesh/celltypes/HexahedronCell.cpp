// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/celltypes/HexahedronCell.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/entities/Vertex.h>

#include <algorithm>
#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------

// UFC: Number of Entities
size_t const HexahedronCell::NE[4][4] = { { 1, 0, 0, 0 },
                                          { 2, 1, 0, 0 },
                                          { 4, 4, 1, 0 },
                                          { 8, 12, 6, 1 } };

// UFC: Vertex Coordinates
real const HexahedronCell::VC[8][3] = { { 0.0, 0.0, 0.0 },
                                        { 1.0, 0.0, 0.0 },
                                        { 1.0, 1.0, 0.0 },
                                        { 0.0, 1.0, 0.0 },
                                        { 0.0, 0.0, 1.0 },
                                        { 1.0, 0.0, 1.0 },
                                        { 1.0, 1.0, 1.0 },
                                        { 0.0, 1.0, 1.0 } };

// UFC: Edge - Incident Vertices
size_t const HexahedronCell::EIV[12][2] = { { 6, 7 },
                                            { 5, 6 },
                                            { 4, 7 },
                                            { 4, 5 },
                                            { 3, 7 },
                                            { 2, 6 },
                                            { 2, 3 },
                                            { 1, 5 },
                                            { 1, 2 },
                                            { 0, 4 },
                                            { 0, 3 },
                                            { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
size_t const HexahedronCell::ENV[12][6] = { { 0, 1, 2, 3, 4, 5 },
                                            { 0, 1, 2, 3, 4, 7 },
                                            { 0, 1, 2, 3, 5, 6 },
                                            { 0, 1, 2, 3, 6, 7 },
                                            { 0, 1, 2, 4, 5, 6 },
                                            { 0, 1, 3, 4, 5, 7 },
                                            { 0, 1, 4, 5, 6, 7 },
                                            { 0, 2, 3, 4, 6, 7 },
                                            { 0, 3, 4, 5, 6, 7 },
                                            { 1, 2, 3, 5, 6, 7 },
                                            { 1, 2, 4, 5, 6, 7 },
                                            { 2, 3, 4, 5, 6, 7 } };

// UFC: Face - Incident Vertices
size_t const HexahedronCell::FIV[6][4] = { { 4, 5, 6, 7 },
                                           { 2, 3, 6, 7 },
                                           { 1, 2, 5, 6 },
                                           { 0, 3, 4, 7 },
                                           { 0, 1, 4, 5 },
                                           { 0, 1, 2, 3 } };

// UFC: Face - Non-Incident Vertices
size_t const HexahedronCell::FNV[6][4] = { { 0, 1, 2, 3 },
                                           { 0, 1, 4, 5 },
                                           { 0, 3, 4, 7 },
                                           { 1, 2, 5, 6 },
                                           { 2, 3, 6, 7 },
                                           { 4, 5, 6, 7 } };

//-----------------------------------------------------------------------------
HexahedronCell::HexahedronCell()
  : CellType( "hexahedron", CellType::hexahedron, CellType::quadrilateral )
{
}
//-----------------------------------------------------------------------------
void HexahedronCell::create_entities( size_t **      e,
                                      size_t         dim,
                                      size_t const * v ) const
{
  // We only need to know how to create edges and faces
  switch ( dim )
  {
    case 1:
      // Create the twelve edges
      e[0][0]  = v[6];
      e[0][1]  = v[7];
      e[1][0]  = v[5];
      e[1][1]  = v[6];
      e[2][0]  = v[4];
      e[2][1]  = v[7];
      e[3][0]  = v[4];
      e[3][1]  = v[5];
      e[4][0]  = v[3];
      e[4][1]  = v[7];
      e[5][0]  = v[2];
      e[5][1]  = v[6];
      e[6][0]  = v[2];
      e[6][1]  = v[3];
      e[7][0]  = v[1];
      e[7][1]  = v[5];
      e[8][0]  = v[1];
      e[8][1]  = v[2];
      e[9][0]  = v[0];
      e[9][1]  = v[4];
      e[10][0] = v[0];
      e[10][1] = v[3];
      e[11][0] = v[0];
      e[11][1] = v[1];
      break;
    case 2:
      // Create the six faces
      e[0][0] = v[4];
      e[0][1] = v[5];
      e[0][2] = v[6];
      e[0][3] = v[7];
      e[1][0] = v[2];
      e[1][1] = v[3];
      e[1][2] = v[7];
      e[1][3] = v[6];
      e[2][0] = v[1];
      e[2][1] = v[2];
      e[2][2] = v[6];
      e[2][3] = v[5];
      e[3][0] = v[0];
      e[3][1] = v[4];
      e[3][2] = v[7];
      e[3][3] = v[3];
      e[4][0] = v[0];
      e[4][1] = v[1];
      e[4][2] = v[5];
      e[4][3] = v[4];
      e[5][0] = v[0];
      e[5][1] = v[3];
      e[5][2] = v[2];
      e[5][3] = v[1];
      break;
    default:
      error( "Invalid topological dimension for creation of entities: %d.",
             dim );
      break;
  }
}
//-----------------------------------------------------------------------------
void HexahedronCell::order_entities( MeshTopology & topology, size_t i ) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1, 3 - 0, 3 - 1, 3 - 2
  dolfin_assert( topology.type().cellType() == this->cell_type );

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if ( topology.connectivity( 1, 0 ) )
  {
    dolfin_assert( topology.connectivity( 3, 1 ) );

    // Get edges
    Array< size_t > const & cell_edges = topology( 3, 1 )[i];
    dolfin_assert( cell_edges.size() >= 12 );

    // Sort vertices on each edge
    for ( size_t i = 0; i < 12; ++i )
    {
      Array< size_t > & edge_vertices = topology( 1, 0 )[cell_edges[i]];
      dolfin_assert( edge_vertices.size() >= 2 );
      std::sort( edge_vertices.data(), edge_vertices.data() + 2 );
    }
  }

  // Sort local vertices on facets in ascending order, connectivity 2 - 0
  if ( topology.connectivity( 2, 0 ) )
  {
    dolfin_assert( topology.connectivity( 3, 2 ) );

    // Get facets
    Array< size_t > const & cell_facets = topology( 3, 2 )[i];
    dolfin_assert( cell_facets.size() >= 6 );

    // Sort vertices on each facet
    for ( size_t i = 0; i < 6; ++i )
    {
      Array< size_t > & facet_vertices = topology( 2, 0 )[cell_facets[i]];
      dolfin_assert( facet_vertices.size() >= 4 );
      std::sort( facet_vertices.data(), facet_vertices.data() + 4 );
    }
  }

  // Sort local edges on local facets after non-incident vertex, connectivity 2
  // - 1
  if ( topology.connectivity( 2, 1 ) )
  {
    dolfin_assert( topology.connectivity( 3, 2 ) );
    dolfin_assert( topology.connectivity( 2, 0 ) );
    dolfin_assert( topology.connectivity( 1, 0 ) );

    // Get facet numbers
    Array< size_t > const & cell_facets = topology( 2, 1 )[i];

    // Loop over facets on cell
    for ( size_t i = 0; i < 6; ++i )
    {
      // For each facet number get the global vertex numbers
      // Facet
      Array< size_t > const & facet_vertices = topology( 2, 0 )[cell_facets[i]];

      // For each facet number get the global edge number
      Array< size_t > & cell_edges = topology( 2, 1 )[cell_facets[i]];

      // Loop over vertices on facet
      size_t m = 0;
      for ( size_t j = 0; j < 4; ++j )
      {
        // Loop edges on facet
        for ( size_t k( m ); k < 4; ++k )
        {
          // For each edge number get the global vertex numbers
          Array< size_t > const & edge_vertices =
            topology( 1, 0 )[cell_edges[k]];

          // Check if the jth vertex of facet i is non-incident on edge k
          if ( !std::count( edge_vertices.data(),
                            edge_vertices.data() + 2,
                            facet_vertices[j] ) )
          {
            // Swap facet numbers
            size_t tmp    = cell_edges[m];
            cell_edges[m] = cell_edges[k];
            cell_edges[k] = tmp;
            m++;
            break;
          }
        }
      }
    }
  }

  // Vertices are supposed to have been ordered appropriately

  // Sort local edges on cell after non-incident vertex tuple, connectivity 3-1
  if ( topology.connectivity( 3, 1 ) )
  {
    dolfin_assert( topology.connectivity( 1, 0 ) );

    // Get cell vertices and edge numbers
    Array< size_t > const & cell_vertices = topology( 3, 0 )[i];
    Array< size_t > &       cell_edges    = topology( 3, 1 )[i];

    // Loop two vertices on cell as a lexicographical tuple
    size_t m = 0;
    for ( size_t t = 0; t < 12; ++t )
    {
      size_t v0 = ENV[t][0];
      size_t v1 = ENV[t][1];
      // Loop edge numbers
      for ( size_t k = m; k < 12; ++k )
      {
        // Get local vertices on edge
        Array< size_t > const & edge_vertices = topology( 1, 0 )[cell_edges[k]];

        // Check if the ith and jth vertex of the cell are non-incident on edge
        // k
        if ( !std::count( edge_vertices.data(),
                          edge_vertices.data() + 2,
                          cell_vertices[v0] )
             && !std::count( edge_vertices.data(),
                             edge_vertices.data() + 2,
                             cell_vertices[v1] ) )
        {
          // Swap edge numbers
          size_t tmp    = cell_edges[m];
          cell_edges[m] = cell_edges[k];
          cell_edges[k] = tmp;
          ++m;
          break;
        }
      }
    }
  }

  // Sort local facets on cell after non-incident vertex, connectivity 3 - 2
  if ( topology.connectivity( 3, 2 ) )
  {
    dolfin_assert( topology.connectivity( 2, 0 ) );

    // Get cell vertices and facet numbers
    Array< size_t > const & cell_vertices = topology( 3, 0 )[i];
    Array< size_t > &       cell_facets   = topology( 3, 2 )[i];

    //
    size_t m = 0;
    for ( size_t t = 0; t < 6; ++t )
    {
      size_t v0 = FNV[t][0];
      size_t v1 = FNV[t][1];
      size_t v2 = FNV[t][2];
      size_t v3 = FNV[t][3];

      // Loop facets on cell
      for ( size_t k = m; k < 6; ++k )
      {
        Array< size_t > const & facet_vertices =
          topology( 2, 0 )[cell_facets[k]];

        // Check if the ith vertex of the cell is non-incident on facet j
        if ( !std::count( facet_vertices.data(),
                          facet_vertices.data() + 4,
                          cell_vertices[v0] )
             && !std::count( facet_vertices.data(),
                             facet_vertices.data() + 4,
                             cell_vertices[v1] )
             && !std::count( facet_vertices.data(),
                             facet_vertices.data() + 4,
                             cell_vertices[v2] )
             && !std::count( facet_vertices.data(),
                             facet_vertices.data() + 4,
                             cell_vertices[v3] ) )
        {
          // Swap facet numbers
          size_t tmp     = cell_facets[m];
          cell_facets[m] = cell_facets[k];
          cell_facets[k] = tmp;
          ++m;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void HexahedronCell::order_facet( size_t vertices[], Facet & facet ) const
{
  static size_t e[4] = { 0 };
  Mesh &        mesh = facet.mesh();
  Cell const    cell( mesh, facet.entities( TD )[0] );

  // Find facet, it is outward oriented by definition.
  // NOTE: Keep consistency with create_entities.
  size_t                  facet_index = cell.index( facet );
  Array< size_t > const & v           = cell.entities( 0 );
  switch ( facet_index )
  {
    case 0:
      e[0] = v[4];
      e[1] = v[5];
      e[2] = v[6];
      e[3] = v[7];
      break;
    case 1:
      e[0] = v[2];
      e[1] = v[3];
      e[2] = v[7];
      e[3] = v[6];
      break;
    case 2:
      e[0] = v[1];
      e[1] = v[2];
      e[2] = v[6];
      e[3] = v[5];
      break;
    case 3:
      e[0] = v[0];
      e[1] = v[4];
      e[2] = v[7];
      e[3] = v[3];
      break;
    case 4:
      e[0] = v[0];
      e[1] = v[1];
      e[2] = v[5];
      e[3] = v[4];
      break;
    case 5:
      e[0] = v[0];
      e[1] = v[3];
      e[2] = v[2];
      e[3] = v[1];
      break;
    default:
      error( "HexahedronCell : invalid local facet index %u", facet_index );
      break;
  }

  // Reorder facet vertices
  for ( size_t i = 0; i < 4; ++i )
  {
    for ( size_t j = 0; j < 4; ++j )
    {
      if ( e[i] == facet.entities( 0 )[j] )
      {
        e[i] = vertices[j];
        break;
      }
    }
  }
  std::memcpy( &vertices[0], &e[0], 4 * sizeof( size_t ) );
}
//-----------------------------------------------------------------------------
void HexahedronCell::refine_cell( Cell &       cell,
                                  MeshEditor & editor,
                                  size_t &     current_cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get vertices and edges
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( !v.empty() );
  Array< size_t > const & e = cell.entities( 1 );
  dolfin_assert( !e.empty() );
  Array< size_t > const & f = cell.entities( 2 );
  dolfin_assert( !f.empty() );

  // Compute indices for the twenty-seven new vertices
  size_t const v00     = v[0];
  size_t const v01     = v[1];
  size_t const v02     = v[2];
  size_t const v03     = v[3];
  size_t const v04     = v[4];
  size_t const v05     = v[5];
  size_t const v06     = v[6];
  size_t const v07     = v[7];
  size_t const eoffset = cell.mesh().size( 0 );
  size_t const e00     = eoffset + e[findEdge( 0, cell )];
  size_t const e01     = eoffset + e[findEdge( 1, cell )];
  size_t const e02     = eoffset + e[findEdge( 2, cell )];
  size_t const e03     = eoffset + e[findEdge( 3, cell )];
  size_t const e04     = eoffset + e[findEdge( 4, cell )];
  size_t const e05     = eoffset + e[findEdge( 5, cell )];
  size_t const e06     = eoffset + e[findEdge( 6, cell )];
  size_t const e07     = eoffset + e[findEdge( 7, cell )];
  size_t const e08     = eoffset + e[findEdge( 8, cell )];
  size_t const e09     = eoffset + e[findEdge( 9, cell )];
  size_t const e10     = eoffset + e[findEdge( 10, cell )];
  size_t const e11     = eoffset + e[findEdge( 11, cell )];
  size_t const foffset = eoffset + cell.mesh().size( 1 );
  size_t const f00     = foffset + f[findFace( 0, cell )];
  size_t const f01     = foffset + f[findFace( 1, cell )];
  size_t const f02     = foffset + f[findFace( 2, cell )];
  size_t const f03     = foffset + f[findFace( 3, cell )];
  size_t const f04     = foffset + f[findFace( 4, cell )];
  size_t const f05     = foffset + f[findFace( 5, cell )];
  size_t const coffset = foffset + cell.mesh().size( 2 );
  size_t const c00     = coffset + cell.index();

  // Add the eight new cells
  size_t const cv0[8] = { v00, e11, f05, e10, e09, f04, c00, f03 };
  editor.add_cell( current_cell++, &cv0[0] );
  size_t const cv1[8] = { e11, v01, e08, f05, f04, e07, f02, c00 };
  editor.add_cell( current_cell++, &cv1[0] );
  size_t const cv2[8] = { f05, e08, v02, e06, c00, f02, e05, f01 };
  editor.add_cell( current_cell++, &cv2[0] );
  size_t const cv3[8] = { e10, f05, e06, v03, f03, c00, f01, e04 };
  editor.add_cell( current_cell++, &cv3[0] );
  size_t const cv4[8] = { e09, f04, c00, f03, v04, e03, f00, e02 };
  editor.add_cell( current_cell++, &cv4[0] );
  size_t const cv5[8] = { f04, e07, f02, c00, e03, v05, e01, f00 };
  editor.add_cell( current_cell++, &cv5[0] );
  size_t const cv6[8] = { c00, f02, e05, f01, f00, e01, v06, e00 };
  editor.add_cell( current_cell++, &cv6[0] );
  size_t const cv7[8] = { f03, c00, f01, e04, e02, f00, e00, v07 };
  editor.add_cell( current_cell++, &cv7[0] );
}
//-----------------------------------------------------------------------------
auto HexahedronCell::intersects( MeshEntity const & e, Point const & ) const
  -> bool
{
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[3][0] );
  MAYBE_UNUSED( e )

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  size_t const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);
   */

  error( "Collision of hexahedron with point not implemented." );

  return true;
}
//-----------------------------------------------------------------------------
auto HexahedronCell::intersects( MeshEntity const & e,
                                 Point const &,
                                 Point const & ) const -> bool
{
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[3][0] );
  MAYBE_UNUSED( e )

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  size_t const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);
   */

  error( "Collision of hexahedron with segment not implemented." );

  return true;
}
//-----------------------------------------------------------------------------
void HexahedronCell::create_reference_cell( Mesh & mesh ) const
{
  MeshEditor me( mesh, CellType::hexahedron, 3, DOLFIN_COMM_SELF );
  me.init_vertices( 8 );
  me.add_vertex( 0, VC[0] );
  me.add_vertex( 1, VC[1] );
  me.add_vertex( 2, VC[2] );
  me.add_vertex( 3, VC[3] );
  me.add_vertex( 4, VC[4] );
  me.add_vertex( 5, VC[5] );
  me.add_vertex( 6, VC[6] );
  me.add_vertex( 7, VC[7] );
  me.init_cells( 1 );
  size_t const cv0[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  me.add_cell( 0, cv0 );
  me.close();
}
//-----------------------------------------------------------------------------
auto HexahedronCell::description() const -> std::string
{
  return std::string( "hexahedron (hypercube of topological dimension 3)" );
}
//-----------------------------------------------------------------------------
void HexahedronCell::disp() const
{
  message( "HexahedronCell" );
  begin( "--------------" );
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
auto HexahedronCell::check( Cell & cell ) const -> bool
{
  bool ret = CellType::check( cell );

  //
  MeshTopology const &    topology = cell.mesh().topology();
  Array< size_t > const & v        = cell.entities( 0 );
  dolfin_assert( !v.empty() );

  // Check edge -> incident vertices mapping
  if ( cell.mesh().topology().connectivity( 1, 0 ) )
  {
    Array< size_t > const & e = cell.entities( 1 );
    dolfin_assert( !e.empty() );
    for ( size_t i = 0; i < 12; ++i )
    {
      Array< size_t > const & ev = topology( 1, 0 )[e[i]];
      dolfin_assert( ev.size() >= 2 );
      for ( size_t j = 0; j < 2; ++j )
      {
        if ( ev[j] != v[EIV[i][0]] && ev[j] != v[EIV[i][1]] )
        {
          ret = false;
          warning( "CellType : invalid edge -> incident vertices mapping\n"
                   "e[%u] = (v%u, v%u)",
                   i,
                   ev[0],
                   ev[1] );
        }
      }
    }
  }

  // Check face -> incident vertices mapping
  if ( cell.mesh().topology().connectivity( 2, 0 ) )
  {
    Array< size_t > const & f = cell.entities( 2 );
    dolfin_assert( f.size() >= 6 );
    for ( size_t i = 0; i < 6; ++i )
    {
      Array< size_t > const & fv = topology( 2, 0 )[f[i]];
      dolfin_assert( fv.size() >= 4 );
      for ( size_t j = 0; j < 4; ++j )
      {
        if ( fv[j] != v[FIV[i][0]] && fv[j] != v[FIV[i][1]]
             && fv[j] != v[FIV[i][2]] && fv[j] != v[FIV[i][3]] )
        {
          ret = false;
          warning( "CellType : invalid face -> incident vertices mapping\n"
                   "f[%u] = (v%u, v%u, v%u, v%u)",
                   i,
                   fv[0],
                   fv[1],
                   fv[2],
                   fv[3] );
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------
auto HexahedronCell::findEdge( size_t i, Cell const & cell ) const -> size_t
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( !v.empty() );
  Array< size_t > const & e = cell.entities( 1 );
  dolfin_assert( !e.empty() );

  // Look for edge satisfying ordering convention
  MeshTopology const & topology = cell.mesh().topology();
  for ( size_t j = 0; j < 12; ++j )
  {
    Array< size_t > const & ev = topology( 1, 0 )[e[j]];
    dolfin_assert( not ev.empty() );
    // Check incident pairs instead of non-incident quadruples
    size_t const v0 = v[EIV[i][0]];
    size_t const v1 = v[EIV[i][1]];
    if ( ( ev[0] == v0 && ev[1] == v1 ) || ( ev[0] == v1 && ev[1] == v0 ) )
    {
      return j;
    }
  }

  // We should not reach this
  error( "Unable to find edge with index %d in hexahedron.", cell.index() );

  return 0;
}
//-----------------------------------------------------------------------------
auto HexahedronCell::findFace( size_t i, Cell const & cell ) const -> size_t
{
  // Ordering convention for faces (order of non-incident vertices)

  // Get vertices and edges
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( !v.empty() );
  Array< size_t > const & f = cell.entities( 2 );
  dolfin_assert( !f.empty() );

  // Look for edge satisfying ordering convention
  MeshTopology const & topology = cell.mesh().topology();
  size_t const         v0       = v[FNV[i][0]];
  size_t const         v1       = v[FNV[i][1]];
  size_t const         v2       = v[FNV[i][2]];
  size_t const         v3       = v[FNV[i][3]];
  for ( size_t j = 0; j < 6; ++j )
  {
    Array< size_t > const & fv = topology( 2, 0 )[f[j]];
    dolfin_assert( not fv.empty() );
    if ( fv[0] != v0 && fv[0] != v1 && fv[0] != v2 && fv[0] != v3 && fv[1] != v0
         && fv[1] != v1 && fv[1] != v2 && fv[1] != v3 && fv[2] != v0
         && fv[2] != v1 && fv[2] != v2 && fv[2] != v3 && fv[3] != v0
         && fv[3] != v1 && fv[3] != v2 && fv[3] != v3 )
    {
      return j;
    }
  }

  // We should not reach this
  error( "Unable to find face with index %d in hexahedron.", cell.index() );

  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
