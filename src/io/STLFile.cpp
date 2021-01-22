// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/STLFile.h>

#include <dolfin/common/byteswap.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/MeshEditor.h>

#include <array>
#include <fstream>

namespace dolfin
{

//--- Implementation detail ---------------------------------------------------
struct stl_vertex
{
  double v[3];
  uint   index;

  auto operator<( const stl_vertex & other ) const -> bool
  {
    return ( ( v[0] < other.v[0] )
             or ( v[0] == other.v[0]
                  and ( v[1] < other.v[1]
                        or ( v[1] == other.v[1] and v[2] < other.v[2] ) ) ) );
  }

  auto operator==( const stl_vertex & other ) const -> bool
  {
    return ( v[0] == other.v[0] && v[1] == other.v[1] && v[2] == other.v[2] );
  }
};

//-----------------------------------------------------------------------------
STLFile::STLFile( const std::string filename )
  : GenericFile( "STL", filename )
  , min( Point( -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX, -DOLFIN_REAL_MAX ) )
  , max( Point( +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX, +DOLFIN_REAL_MAX ) )
{
}
//-----------------------------------------------------------------------------
void STLFile::set_min( Point const & min_ )
{
  min = min_;
}
//-----------------------------------------------------------------------------
void STLFile::set_max( Point const & max_ )
{
  max = max_;
}
//-----------------------------------------------------------------------------
void STLFile::operator>>( Mesh & mesh )
{
  std::array< char, 80 > hdr;
  std::array< float, 3 > data;
  uint  ntri = 0;

  struct stl_vertex              V[3];
  _ordered_set< stl_vertex >     vertices;
  std::vector< std::array< size_t, 3 > > cells;

  std::ifstream fp( filename.c_str(), std::ifstream::binary );
  fp.read( hdr.data(), hdr.size() * sizeof( char ) );
  fp.read( reinterpret_cast< char * >( &ntri ), sizeof( uint ) );

#ifdef HAVE_BIG_ENDIAN
  ntri = bswap( ntri );
#endif

  uint v_index = 0;
  uint c_index = 0;
  for ( uint i = 0; i < ntri; ++i )
  {
    /* Normal */
    fp.read( reinterpret_cast< char * >( &data ), 3 * sizeof( float ) );

    // load data from file
    for ( uint j = 0; j < 3; ++j )
    {
      /* Vertex v1 v2 v3 */
      fp.read( reinterpret_cast< char * >( &data ), 3 * sizeof( float ) );

#ifdef HAVE_BIG_ENDIAN
      V[j].v[0] = static_cast< double >( bswap( data[0] ) );
      V[j].v[1] = static_cast< double >( bswap( data[1] ) );
      V[j].v[2] = static_cast< double >( bswap( data[2] ) );
#else
      V[j].v[0] = static_cast< double >( data[0] );
      V[j].v[1] = static_cast< double >( data[1] );
      V[j].v[2] = static_cast< double >( data[2] );
#endif
    }

    // check if at least one of the three vertices are inside the bounding box
    bool inside = false;

    for ( uint j = 0; j < 3; ++j )
    {
      if (     ( V[j].v[0] >= min[0] and V[j].v[0] <= max[0] )
           and ( V[j].v[1] >= min[1] and V[j].v[1] <= max[1] )
           and ( V[j].v[2] >= min[2] and V[j].v[2] <= max[2] ) )
      {
        inside = true;
      }
    }

    // add vertices and cells to the list
    if ( inside == true )
    {
      cells.resize( cells.size() + 1 );

      for ( uint j = 0; j < 3; ++j )
      {
        if ( vertices.find( V[j] ) != vertices.end() )
        {
          cells.back()[j] = vertices.find( V[j] )->index;
        }
        else
        {
          V[j].index      = v_index++;
          cells.back()[j] = V[j].index;
          vertices.insert( V[j] );
        }
      }
    }

    // Aux data
    fp.read( reinterpret_cast< char * >( &hdr ), 2 * sizeof( char ) );
  }

  fp.close();

  MeshEditor editor( mesh, CellType::triangle, 3, MPI::DOLFIN_COMM_SELF );

  // add vertices
  editor.init_cells( cells.size() );
  for ( uint i = 0; i < cells.size(); ++i )
  {
    editor.add_cell( c_index++, &cells[i][0] );
  }

  // add cells
  editor.init_vertices( vertices.size() );
  for ( stl_vertex const & vert : vertices )
  {
    editor.add_vertex( vert.index, &vert.v[0] );
  }

  editor.close();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
