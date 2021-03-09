// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/VTKFile.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/Encoder.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

#include <cstdint>

namespace dolfin
{

//----------------------------------------------------------------------------

VTKFile::VTKFile( const std::string filename )
  : GenericFile( "VTK", filename )
  , _t( nullptr )
{
}

//----------------------------------------------------------------------------

VTKFile::VTKFile( const std::string filename, real const & t )
  : GenericFile( "VTK", filename )
  , _t( &t )
{
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( Mesh & mesh )
{
  // Update vtu file name and clear file
  vtuNameUpdate( counter );

  // Only the root updates the pvd file
  if ( mesh.is_distributed() )
  {
    if ( PE::rank() == 0 )
    {
      // Update pvtu file name and clear file
      pvtuNameUpdate( counter );

      // Write pvtu file
      pvtuFileWrite( false, 0 );

      // Write pvd file
      pvdFileWrite( counter, true );
    }
  }
  else
  {
#if DOLFIN_HAVE_MPI
    if ( PE::rank() == 0 && MPI::size() > 1 )
    {
      warning( "Writing serial mesh in a parallel run" );
    }
#endif

    // Write pvd file
    pvdFileWrite( counter, false );
  }

  // Write headers
  VTKHeaderOpen( mesh );

  // Write mesh
  MeshWrite( mesh );

  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the mesh
  ++counter;

  message( 1,
           "VTKFile: Saved mesh %s (%s) to file %s in VTK format.",
           mesh.name().c_str(),
           mesh.label().c_str(),
           filename.c_str() );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( MeshFunction< int > & meshfunction )
{
  MeshFunctionWrite( meshfunction );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( MeshFunction< size_t > & meshfunction )
{
  MeshFunctionWrite( meshfunction );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( MeshFunction< real > & meshfunction )
{
  MeshFunctionWrite( meshfunction );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( MeshFunction< bool > & meshfunction )
{
  MeshFunctionWrite( meshfunction );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( Function & u )
{
  LabelList< Function > tmp( 1, Label< Function >( u, "U" ) );
  write_dataset( tmp );
}

//----------------------------------------------------------------------------

void VTKFile::operator<<( LabelList< Function > & list )
{
  write_dataset( list );
}

//----------------------------------------------------------------------------

void VTKFile::read()
{
  opened_read = true;
}

//----------------------------------------------------------------------------

void VTKFile::write()
{
  if ( !opened_write )
  {
    if ( PE::rank() == 0 )
    {
      // Clear file
      FILE * fp = fopen( filename.c_str(), "w" );
      fclose( fp );
    }
  }
  opened_write = true;
}

//----------------------------------------------------------------------------

void VTKFile::write_dataset( LabelList< Function > & f )
{
  dolfin_assert( f.size() > 0 );

  // Update vtu file name and clear file
  vtuNameUpdate( counter );

  /// @todo This only writes the first mesh encountered
  Mesh & mesh = f[0].first->mesh();

  // Only the root updates the pvd file
  if ( mesh.is_distributed() )
  {
    if ( PE::rank() == 0 )
    {
      // Update pvtu file name and clear file
      pvtuNameUpdate( counter );

      // Write pvtu file
      pvtuFileWriteFunction( f );

      // Write pvd file
      pvdFileWrite( counter, true );
    }
  }
  else
  {
#if DOLFIN_HAVE_MPI
    if ( PE::rank() == 0 && MPI::size() > 1 )
    {
      warning( "Writing serial mesh in a parallel run" );
    }
#endif

    // Write pvd file
    pvdFileWrite( counter, false );
  }

  // Write headers
  VTKHeaderOpen( mesh );

  // Write Mesh
  MeshWrite( mesh );

  // Write results
  ResultsWrite( f );

  // Close headers
  VTKHeaderClose();

  // Increase the number of times we have saved the function
  ++counter;

  message( 1,
           "VTKFile: Saved %d %s to file %s in VTK format.",
           f.size(),
           ( f.size() > 1 ? "functions" : "function" ),
           filename.c_str() );
}

//----------------------------------------------------------------------------

void VTKFile::MeshWrite( Mesh & mesh ) const
{
  // Open file
  FILE * fp = fopen( vtu_filename.c_str(), "a" );

  size_t const d                     = mesh.geometry_dimension();
  size_t const num_mesh_verts        = mesh.size( 0 );
  size_t const num_mesh_cells        = mesh.num_cells();
  size_t const num_cell_verts        = mesh.type().num_entities( 0 );
  size_t const cell_verts_block_size = num_cell_verts * num_mesh_cells;

  // Write vertex positions
  fprintf( fp, "<Points>\n" );
  fprintf( fp,
           "<DataArray  type=\"Float32\"  NumberOfComponents=\"3\"  "
           "format=\"binary\">\n" );
  std::vector< float >           v_data( 3 * num_mesh_verts, 0.0 );
  std::vector< float >::iterator entry = v_data.begin();
  for ( VertexIterator v( mesh ); !v.end(); ++v, entry += 3 )
  {
    std::copy( v->x(), v->x() + d, entry );
  }

  // Create encoded stream
  std::stringstream base64_stream;
  encode_stream( base64_stream, v_data );
  fprintf( fp, "%s\n", base64_stream.str().c_str() );
  fprintf( fp, "</DataArray>\n" );
  fprintf( fp, "</Points>\n" );

  // Write cell connectivity
  fprintf( fp, "<Cells>\n" );
  fprintf(
    fp,
    "<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"binary\">\n" );
  std::vector< int32_t >           c_data( cell_verts_block_size );
  std::vector< int32_t >::iterator c_entry = c_data.begin();
  for ( CellIterator c( mesh ); !c.end(); ++c )
  {
    for ( VertexIterator v( *c ); !v.end(); ++v )
    {
      *c_entry++ = static_cast< int32_t >( v->index() );
    }
  }

  // Create encoded stream
  std::stringstream base64_c_stream;
  encode_stream( base64_c_stream, c_data );
  fprintf( fp, "%s\n", base64_c_stream.str().c_str() );
  fprintf( fp, "</DataArray>\n" );

  // Write offset into connectivity array for the end of each cell
  fprintf(
    fp, "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"binary\">\n" );
  std::vector< int32_t >           o_data( cell_verts_block_size );
  std::vector< int32_t >::iterator oo_entry = o_data.begin();
  for ( size_t offsets = 1; offsets <= num_mesh_cells; ++offsets )
  {
    *oo_entry++ = static_cast< int32_t >( offsets * num_cell_verts );
  }

  std::stringstream base64_oo_stream;
  encode_stream( base64_oo_stream, o_data );
  fprintf( fp, "%s\n", base64_oo_stream.str().c_str() );
  fprintf( fp, "</DataArray>\n" );

  // Write cell type
  fprintf( fp,
           "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"binary\">\n" );

  uint8_t typeval = 0;
  switch ( mesh.type().cellType() )
  {
    case CellType::hexahedron:
      typeval = uint8_t( 12 );
      break;
    case CellType::tetrahedron:
      typeval = uint8_t( 10 );
      break;
    case CellType::quadrilateral:
      typeval = uint8_t( 9 );
      break;
    case CellType::triangle:
      typeval = uint8_t( 5 );
      break;
    case CellType::interval:
      typeval = uint8_t( 3 );
      break;
    case CellType::point:
      typeval = uint8_t( 1 );
      break;
    default:
      error( "Provided mesh type is not supported" );
      break;
  }
  std::vector< uint8_t > t_data( num_mesh_cells, typeval );

  // Create encoded stream
  std::stringstream base64_t_stream;
  encode_stream( base64_t_stream, t_data );
  fprintf( fp, "%s\n", base64_t_stream.str().c_str() );
  fprintf( fp, "</DataArray>\n" );
  fprintf( fp, "</Cells>\n" );

  // Close file
  fclose( fp );
}

//----------------------------------------------------------------------------

void VTKFile::ResultsWrite( LabelList< Function > f ) const
{
  // Open file
  FILE * fp = fopen( vtu_filename.c_str(), "a" );

  // Assume same mesh for all data arrays
  Mesh & mesh = f[0].first->mesh();
  for ( LabelList< Function >::iterator it = f.begin(); it != f.end(); it++ )
  {
    if ( it->first->mesh() != mesh )
    {
      error( "VTKFile implementation supports only one mesh per file." );
    }
  }

  //--- Interpolation to vertices for general functions-----------------------
  fprintf( fp, "<PointData>\n" );
  for ( LabelList< Function >::iterator it = f.begin(); it != f.end(); it++ )
  {
    Function * u = it->first;

    // Check type of function space
    if ( u->space().is_cellwise_constant() )
    {
      // These are cell based function and will be written as CellData
      continue;
    }

    //
    size_t const value_rank = u->rank();
    size_t const value_dim  = u->dim( 0 );
    size_t const value_size = u->value_size();

    // Write function data at mesh vertices
    std::string & name = it->second;
    switch ( value_rank )
    {
      case 0:
        fprintf(
          fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  format=\"binary\">\n",
              name.c_str() );
        break;
      case 1:
        fprintf( fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  "
                     "NumberOfComponents=\"%ld\" format=\"binary\">\n",
                     name.c_str(), std::max( value_dim, ( size_t ) 3 ) );
        break;
      default:
        fprintf( fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  "
                     "NumberOfComponents=\"%ld\" format=\"binary\">\n",
                     name.c_str(), value_size );
        break;
    }

    // Allocate memory for function values at vertices
    size_t const num_verts = mesh.size( 0 );
    real *       values    = new real[value_size * num_verts];

    // Get function values at vertices
    u->interpolate_vertex_values( values );

    // Copy values
    // Should be value_dim^value_rank but 1^0 = 1 and (2|3|n)^1 = (2|3|n) so ...
    std::vector< float > data;

    switch ( value_size )
    {
      case 1:
        data.resize( num_verts );
        {
          std::vector< float >::iterator entry = data.begin();
          for ( VertexIterator vertex( mesh ); !vertex.end(); ++vertex )
          {
            *entry++ = values[vertex->index()];
          }
        }
        break;
      case 2:
        data.resize( 3 * num_verts );
        {
          std::vector< float >::iterator entry = data.begin();
          for ( VertexIterator vertex( mesh ); !vertex.end(); ++vertex )
          {
            *entry++ = values[vertex->index()];
            *entry++ = values[vertex->index() + num_verts];
            *entry++ = 0.0;
          }
        }
        break;
      case 3:
        data.resize( 3 * num_verts );
        {
          std::vector< float >::iterator entry = data.begin();
          for ( VertexIterator vertex( mesh ); !vertex.end(); ++vertex )
          {
            *entry++ = values[vertex->index()];
            *entry++ = values[vertex->index() + num_verts];
            *entry++ = values[vertex->index() + 2 * num_verts];
          }
        }
        break;
      default:
        data.resize( value_size * num_verts );
        {
          std::vector< float >::iterator entry = data.begin();
          for ( VertexIterator vertex( mesh ); !vertex.end(); ++vertex )
          {
            for ( size_t d = 0; d < value_size; ++d )
            {
              *entry++ = values[vertex->index() + d * num_verts];
            }
          }
        }
        break;
    }

    //
    delete[] values;

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream( base64_stream, data );
    fprintf( fp, "%s\n", base64_stream.str().c_str() );
    fprintf( fp, "</DataArray>\n" );
  }
  fprintf( fp, "</PointData>\n" );

  //--- Cellwise functions----------------------------------------------------
  fprintf( fp, "<CellData>\n" );
  for ( LabelList< Function >::iterator it = f.begin(); it != f.end(); it++ )
  {
    Function * u = it->first;

    // Check type of function space
    if ( !u->space().is_cellwise_constant() )
    {
      // These are not cell based functions and will be written as PointData
      continue;
    }

    //
    size_t const value_rank = u->rank();
    size_t const value_dim  = u->dim( 0 );
    size_t const value_size = u->value_size();

    // Write function data in cell
    std::string &        name = it->second;
    std::vector< float > data;
    switch ( value_rank )
    {
      case 0:
        fprintf( fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  format=\"binary\">\n",
                     name.c_str() );
        break;
      case 1:
        fprintf( fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  "
                     "NumberOfComponents=\"%ld\" format=\"binary\">\n",
                     name.c_str(), std::max( value_dim, ( size_t ) 3 ) );
        break;
      default:
        fprintf( fp, "<DataArray  type=\"Float32\"  Name=\"%s\"  "
                     "NumberOfComponents=\"%ld\" format=\"binary\">\n",
                     name.c_str(), value_size );
        break;
    }

    // Copy values
    // Should be value_dim^value_rank but 1^0 = 1 and (2|3|n)^1 = (2|3|n) so ...
    real * values = nullptr;
    u->get_block( values );
    switch ( value_size )
    {
      case 1:
        data.resize( mesh.num_cells() );
        {
          std::vector< float >::iterator entry = data.begin();
          size_t                         ii    = 0;
          for ( CellIterator cell( mesh ); !cell.end(); ++cell, ++ii )
          {
            *entry++ = values[ii];
          }
        }
        break;
      case 2:
        data.resize( 3 * mesh.num_cells() );
        {
          std::vector< float >::iterator entry = data.begin();
          size_t                         ii    = 0;
          for ( CellIterator cell( mesh ); !cell.end(); ++cell )
          {
            *entry++ = values[ii++];
            *entry++ = values[ii++];
            *entry++ = 0.0;
          }
        }
        break;
      case 3:
        data.resize( 3 * mesh.num_cells() );
        {
          std::vector< float >::iterator entry = data.begin();
          size_t                         ii    = 0;
          for ( CellIterator cell( mesh ); !cell.end(); ++cell )
          {
            *entry++ = values[ii++];
            *entry++ = values[ii++];
            *entry++ = values[ii++];
          }
        }
        break;
      default:
        data.resize( value_size * mesh.num_cells() );
        {
          std::vector< float >::iterator entry = data.begin();
          size_t                         ii    = 0;
          for ( VertexIterator vertex( mesh ); !vertex.end(); ++vertex )
          {
            for ( size_t d = 0; d < value_size; ++d )
            {
              *entry++ = values[ii++];
            }
          }
        }
        break;
    }

    //
    delete[] values;

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream( base64_stream, data );
    fprintf( fp, "%s\n", base64_stream.str().c_str() );
    fprintf( fp, "</DataArray>\n" );
  }
  fprintf( fp, "</CellData>\n" );

  // Close file
  fclose( fp );
}

//----------------------------------------------------------------------------

void VTKFile::pvdFileWrite( size_t num, bool parallel )
{
  std::fstream pvdFile;

  if ( num == 0 )
  {
    // Open pvd file
    pvdFile.open( filename.c_str(), std::ios::out | std::ios::trunc );
    // Write header
    pvdFile << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" >\n"
            << "<Collection>\n";
  }
  else
  {
    // Open pvd file
    pvdFile.open( filename.c_str(), std::ios::out | std::ios::in );
    pvdFile.seekp( mark );
  }

  // Remove directory path from name for pvd file
  std::string fname;
#if DOLFIN_HAVE_MPI
  if ( parallel )
  {
    fname.assign(
      pvtu_filename, filename.find_last_of( "/" ) + 1, pvtu_filename.size() );
  }
  else
  {
    fname.assign(
      vtu_filename, filename.find_last_of( "/" ) + 1, vtu_filename.size() );
  }
#else
  MAYBE_UNUSED( parallel )
  fname.assign(
    vtu_filename, filename.find_last_of( "/" ) + 1, vtu_filename.size() );
#endif

  // Data file name
  if ( _t )
  {
    pvdFile << "<DataSet timestep=\"" << *_t << "\" part=\"0\""
            << " file=\"" << fname << "\"/>\n";
  }
  else
  {
    pvdFile << "<DataSet timestep=\"" << num << "\" part=\"0\""
            << " file=\"" << fname << "\"/>\n";
  }
  mark = pvdFile.tellp();

  // Close headers
  pvdFile << "</Collection>\n";
  pvdFile << "</VTKFile>\n";

  // Close file
  pvdFile.close();
}

//----------------------------------------------------------------------------

void VTKFile::pvtuFileWrite( bool mesh_function, size_t const dim )
{
  std::fstream pvtuFile;

  // Open pvtu file
  pvtuFile.open( pvtu_filename.c_str(), std::ios::out | std::ios::trunc );
  // Write header
  pvtuFile << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n"
           << "<PUnstructuredGrid GhostLevel=\"0\">\n"
           << "<PCellData>\n"
           << "<PDataArray  type=\"Int32\"  Name=\"connectivity\" />\n"
           << "<PDataArray  type=\"Int32\"  Name=\"offsets\" />\n"
           << "<PDataArray  type=\"UInt8\"  Name=\"types\" />\n"
           << "</PCellData>\n"
           << "<PPoints>\n"
           << "<PDataArray  type=\"Float32\"  NumberOfComponents=\"3\" />\n"
           << "</PPoints>\n";

  if ( mesh_function )
  {
    if ( dim > 0 )
    {
      // Cell-based mesh function
      pvtuFile << "<PCellData Scalars=\"U\">\n"
               << "<PDataArray  type=\"Float32\"  Name=\"U\" />\n"
               << "</PCellData>\n";
    }
    else
    {
      // Vertex-based mesh function
      pvtuFile << "<PPointData Scalars=\"U\">\n"
               << "<PDataArray  type=\"Float32\"  Name=\"U\" />\n"
               << "</PPointData>\n";
    }
  }

  std::string fname;
  // Remove rank from vtu filename ( <rank>.vtu)
  fname.assign( vtu_filename, 0, vtu_filename.rfind( "_" ) + 1 );
  if ( fname.rfind( "/" ) != std::string::npos )
    fname = fname.substr( fname.rfind( "/" ) + 1, std::string::npos );

  for ( size_t i = 0; i < MPI::size(); ++i )
    pvtuFile << "<Piece Source=\"" << fname << i << ".vtu\"/>\n";

  pvtuFile << "</PUnstructuredGrid>\n"
           << "</VTKFile>\n";
  pvtuFile.close();

}

//----------------------------------------------------------------------------

void VTKFile::pvtuFileWriteFunction( LabelList< Function > f )
{
  std::fstream pvtuFile;

  // Open pvtu file
  pvtuFile.open( pvtu_filename.c_str(), std::ios::out | std::ios::trunc );
  // Write header
  pvtuFile << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n"
           << "<PUnstructuredGrid GhostLevel=\"0\">\n"
           << "<PPointData>\n";
  for ( LabelList< Function >::iterator it = f.begin(); it != f.end(); it++ )
  {
    Function *    u    = it->first;
    std::string & name = it->second;

    // Check type of function space
    if ( u->space().is_cellwise_constant() )
    {
      // These are cell based function and will be written as CellData
      continue;
    }

    if ( u->rank() == 0 )
    {
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\" />\n";
    }
    else
    {
      // Get number of components
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name
               << "\"  NumberOfComponents=\"" << std::max( u->dim( 0 ), 3ul )
               << "\" />\n";
    }
  }
  pvtuFile << "</PPointData>\n";

  //
  pvtuFile << "<PCellData>\n"
           << "<PDataArray  type=\"Int32\"  Name=\"connectivity\"  />\n"
           << "<PDataArray  type=\"Int32\"  Name=\"offsets\" />\n"
           << "<PDataArray  type=\"UInt8\"  Name=\"types\" />\n"
           << "</PCellData>\n";

  //
  pvtuFile << "<PPoints>\n"
           << "<PDataArray  type=\"Float32\"  NumberOfComponents=\"3\" />\n"
           << "</PPoints>\n";

  //
  pvtuFile << "<PCellData>\n";
  for ( LabelList< Function >::iterator it = f.begin(); it != f.end(); it++ )
  {
    Function *    u    = it->first;
    std::string & name = it->second;

    // Check type of function space
    if ( !u->space().is_cellwise_constant() )
    {
      // These are not cell based functions and will be written as PointData
      continue;
    }

    if ( u->rank() == 0 )
    {
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name << "\" />\n";
    }
    else
    {
      // Get number of components
      pvtuFile << "<PDataArray  type=\"Float32\"  Name=\"" << name
               << "\"  NumberOfComponents=\"" << std::max( u->dim( 0 ), 3ul )
               << "\" />\n";
    }
  }
  pvtuFile << "</PCellData>\n";

  std::string fname;

  // Remove rank from vtu filename ( <rank>.vtu)
  fname.assign( vtu_filename, 0, vtu_filename.rfind( "_" ) + 1 );
  if ( fname.rfind( "/" ) != std::string::npos )
    fname = fname.substr( fname.rfind( "/" ) + 1, std::string::npos );

  for ( size_t i = 0; i < MPI::size(); ++i )
    pvtuFile << "<Piece Source=\"" << fname << i << ".vtu\"/>\n";

  pvtuFile << "</PUnstructuredGrid>\n"
           << "</VTKFile>\n";
  pvtuFile.close();
}

//----------------------------------------------------------------------------

void VTKFile::VTKHeaderOpen( Mesh & mesh ) const
{
  // Open file
  FILE * fp = fopen( vtu_filename.c_str(), "a" );

  // Figure out endianness of machine
  std::string endianness = "";
#if defined HAVE_BIG_ENDIAN
  endianness = "byte_order=\"BigEndian\"";
  ;
#else
  endianness = "byte_order=\"LittleEndian\"";
#endif

  std::string compressor = "";
#ifdef HAVE_LIBZ
  compressor = "compressor=\"vtkZLibDataCompressor\"";
#endif

  // Write headers
  fprintf( fp,
           "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\" %s %s>\n",
           endianness.c_str(),
           compressor.c_str() );
  fprintf( fp, "<UnstructuredGrid>\n" );
  fprintf( fp,
           "<Piece  NumberOfPoints=\" %8lu\"  NumberOfCells=\" %8lu\">\n",
           mesh.size( 0 ),
           mesh.num_cells() );

  // Close file
  fclose( fp );
}

//----------------------------------------------------------------------------

void VTKFile::VTKHeaderClose() const
{
  // Open file
  FILE * fp = fopen( vtu_filename.c_str(), "a" );

  // Close headers
  fprintf( fp, "</Piece>\n </UnstructuredGrid>\n </VTKFile>" );

  // Close file
  fclose( fp );
}

//----------------------------------------------------------------------------

void VTKFile::vtuNameUpdate( const int counter )
{
  std::string        filestart;
  std::ostringstream fileid, newfilename;

  fileid.fill( '0' );
  fileid.width( 6 );

  filestart.assign( filename, 0, filename.rfind( "." ) );

  fileid << counter;
  newfilename << filestart << fileid.str() << "_" << PE::rank() << ".vtu";
  vtu_filename = newfilename.str();

  // Make sure file is empty
  FILE * fp = fopen( vtu_filename.c_str(), "w" );
  fclose( fp );
}

//----------------------------------------------------------------------------

void VTKFile::pvtuNameUpdate( const int counter )
{
  std::string        filestart;
  std::ostringstream fileid, newfilename;

  fileid.fill( '0' );
  fileid.width( 6 );

  filestart.assign( filename, 0, filename.rfind( "." ) );

  fileid << counter;
  newfilename << filestart << fileid.str() << ".pvtu";
  pvtu_filename = newfilename.str();

  // Make sure file is empty
  FILE * fp = fopen( pvtu_filename.c_str(), "w" );
  fclose( fp );
}

//----------------------------------------------------------------------------

template < class T >
void VTKFile::MeshFunctionWrite( MeshFunction< T > & meshfunction )
{
  // Update vtu file name and clear file
  vtuNameUpdate( counter );

  // Write pvd file
  if ( PE::rank() == 0 )
  {
    pvtuNameUpdate( counter );
    pvdFileWrite( counter, meshfunction.mesh().is_distributed() );
    pvtuFileWrite( true, meshfunction.dim() );
  }

  Mesh & mesh = meshfunction.mesh();

  if ( meshfunction.dim() == mesh.topology_dimension() )
  {
    // Write headers
    VTKHeaderOpen( mesh );

    // Write mesh
    MeshWrite( mesh );

    std::vector< float >           data( mesh.num_cells() );
    std::vector< float >::iterator entry = data.begin();
    for ( CellIterator cell( mesh ); !cell.end(); ++cell )
    {
      *entry++ = static_cast< float >( meshfunction( cell->index() ) );
    }

    // Open file
    FILE * fp = fopen( vtu_filename.c_str(), "a" );
    fprintf( fp, "<CellData  Scalars=\"U\">\n" );
    fprintf( fp,
             "<DataArray  type=\"Float32\"  Name=\"U\"  format=\"binary\">\n" );

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream( base64_stream, data );
    fprintf( fp, "%s\n", base64_stream.str().c_str() );
    fprintf( fp, "</DataArray>\n" );
    fprintf( fp, "</CellData>\n" );

    // Close file
    fclose( fp );

    // Close headers
    VTKHeaderClose();
  }
  else if ( meshfunction.dim() == 0 )
  {
    // Write headers
    VTKHeaderOpen( mesh );

    // Write mesh
    MeshWrite( mesh );

    std::vector< float > data;
    data.resize( mesh.size( 0 ) );
    std::vector< float >::iterator entry = data.begin();

    for ( VertexIterator v( mesh ); !v.end(); ++v )
    {
      *entry++ = static_cast< float >( meshfunction( v->index() ) );
    }

    // Open file
    FILE * fp = fopen( vtu_filename.c_str(), "a" );
    fprintf( fp, "<PointData  Scalars=\"U\">\n" );
    fprintf( fp,
             "<DataArray  type=\"Float32\"  Name=\"U\"  format=\"binary\">\n" );

    // Create encoded stream
    std::stringstream base64_stream;
    encode_stream( base64_stream, data );
    fprintf( fp, "%s\n", base64_stream.str().c_str() );
    fprintf( fp, "</DataArray>\n" );
    fprintf( fp, "</PointData>\n" );

    // Close file
    fclose( fp );

    // Close headers
    VTKHeaderClose();
  }
  else
  {
    error( "VTK output of mesh functions is implemented for cell-based and "
           "vertex-based functions only." );
  }

  // Increase the number of times we have saved the mesh function
  ++counter;

  message( 1, "VTKFile: Saved mesh function %d times.", counter );
  message( 1,
           "VTKFile: Saved mesh function ( %s ) to file %s in VTK format.",
           mesh.label().c_str(),
           filename.c_str() );
}

//-----------------------------------------------------------------------------

template < typename T >
void VTKFile::encode_stream( std::stringstream &      stream,
                             const std::vector< T > & data ) const
{

#ifdef HAVE_LIBZ
  encode_inline_compressed_base64( stream, data );
#else
  warning( "zlib must be configured to enable compressed VTK output. "
           "Using uncompressed base64 encoding instead." );
  encode_inline_base64( stream, data );
#endif
}

//----------------------------------------------------------------------------

template < typename T >
void VTKFile::encode_inline_base64( std::stringstream &      stream,
                                    const std::vector< T > & data ) const
{
  const uint32_t size = data.size() * sizeof( T );
  Encoder::encode_base64( &size, 1, stream );
  Encoder::encode_base64( data, stream );
}

//----------------------------------------------------------------------------

#ifdef HAVE_LIBZ
template < typename T >
void VTKFile::encode_inline_compressed_base64(
  std::stringstream &      stream,
  const std::vector< T > & data ) const
{
  uint32_t header[4];
  header[0] = 1;
  header[1] = data.size() * sizeof( T );
  header[2] = 0;

  // Compress data
  std::pair< unsigned char *, unsigned long > compressed_data =
    Encoder::compress_data( data );

  // Length of compressed data
  header[3] = compressed_data.second;

  // Encode header
  Encoder::encode_base64( &header[0], 4, stream );

  // Encode data
  Encoder::encode_base64(
    compressed_data.first, compressed_data.second, stream );

  // Free compress data
  delete[] compressed_data.first;
}
#endif

//----------------------------------------------------------------------------

} // namespace dolfin
