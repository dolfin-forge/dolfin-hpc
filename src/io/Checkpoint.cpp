// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/Checkpoint.h>

#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/ParameterSystem.h>

#include <fstream>
#include <sstream>
#include <string>

namespace dolfin
{

/* Checkpoint file format:
 * 1. Header
 * 1a. time (float)
 * 1b. offset_parametersystem (uint)
 * 1c. offset_mesh (uint)
 * 1d. offset_functions (uint)
 * 1e. offset_vectors (uint)
 * 2. ParameterSystem
 * 3. Mesh
 * 3a. MeshHeader (struct)
 * 3b. Mesh
 * 4. Functions
 * 4a. FunctionsHeader (struct)
 * 4b. Functions
 * 5. Vectors
 * 5a. VectorsHeader (struct)
 * 5b. Vectors
 */

// TODO seekp() for non MPIIO functions

//-----------------------------------------------------------------------------

Checkpoint::Checkpoint()
  : n_( 0 )
  , byte_offset( 0 )
{
}

//-----------------------------------------------------------------------------

Checkpoint::~Checkpoint()
{
}

//-----------------------------------------------------------------------------

void Checkpoint::write( std::string filename, real const t, Mesh & mesh,
                        LabelList< Function > & func,
                        LabelList< GenericVector > & vec )
{
  message( "Writing checkpoint (%s %d) at time %g",
           filename.c_str(), n_, t );

  // deserialize ParameterSystem
  std::string parameters  = ParameterSystem::parameters.serialize();
  uint        param_size  = parameters.size() * sizeof( char );
  chkp_header.offset_mesh = param_size;

  MeshHeader      mesh_header;
  FunctionsHeader functions_header;
  VectorsHeader   vectors_header;

  fill_headers( t, mesh, mesh_header,
                func, functions_header, vec, vectors_header );

#ifndef ENABLE_MPIIO
  stream_t file = stream_t( get_filename( filename ), std::ofstream::binary );

  // write header
  file.write( static_cast< char * >( &chkp_header ), sizeof( CheckpointHeader ) );

  // write ParameterSystem
  file.write( static_cast< char * >( &param_size ), sizeof( uint ) );
  file.write( static_cast< char * >( parameters.c_str() ), param_size );
#else
  stream_t file = MPI_File();
  MPI::file_open( file, get_filename( filename ),
                  MPI_MODE_WRONLY | MPI_MODE_CREATE );

  // write header
  byte_offset = MPI::file_write_all( file, chkp_header, sizeof( CheckpointHeader ) );

  // write ParameterSystem
  byte_offset += MPI::file_write_all( file, param_size );
  byte_offset += MPI::file_write_all( file, parameters.c_str()[0], param_size );
#endif

  write_mesh( file, byte_offset, mesh, mesh_header );
  write_func( file, byte_offset, func, functions_header );
  write_vect( file, byte_offset, vec, vectors_header );

  // increment checkpoint id
  n_++;

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load_parametersystem( std::string filename )
{
  stream_t file = open_file( filename );

  uint param_size = 0;

#ifdef ENABLE_MPIIO
  // load ParameterSystem
  byte_offset += MPI::file_read_all( file, param_size );
  Array< char > p( param_size / sizeof( char ) );
  byte_offset += MPI::file_read_all( file, p[0], param_size );
#else
  // load ParameterSystem
  file.read( static_cast< char * >( &param_size ), sizeof( uint ) );
  Array< char > p( param_size / sizeof( char ) );
  file.read( static_cast< char * >( p.data() ), param_size );
#endif

  message( "Restarting from time %g", chkp_header.time );
  ParameterSystem::parameters.deserialize( std::string( p.data(), p.size() ) );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, Mesh & mesh )
{
  stream_t file = open_file( filename );

  Checkpoint::offset_t byte_offset = chkp_header.offset_mesh;

  MeshHeader mesh_header;

  // load header
#ifdef ENABLE_MPIIO
  MPI::file_read_at_all( file, mesh_header,
                     byte_offset + MPI::rank() * sizeof( MeshHeader ) );
  byte_offset += MPI::size() * sizeof( MeshHeader );
#else
  file.read( ( char * ) &mesh_header, sizeof( MeshHeader ) );
#endif

  Mesh _mesh;
  MeshEditor editor( _mesh, mesh_header.type, mesh_header.gdim );
  editor.init_vertices( mesh_header.num_vertices );

  // load coords
  {
    Array< real > coords( mesh_header.num_coords, 0.0 );

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all(
                     file, coords.data(), mesh_header.num_coords,
                     byte_offset + mesh_header.offsets[0] * sizeof( real ),
                     mesh_header.disp[0] );
 #else
    file.read( ( char * ) coords, ( mesh_header.num_coords ) * sizeof( real ) );
 #endif

    for ( uint i = 0; i < mesh_header.num_coords; i += mesh_header.gdim )
    {
      editor.add_vertex( i / mesh_header.gdim, coords.data() + i );
    }
  }

  // load cells
  {
    editor.init_cells( mesh_header.num_cells );
    Array< uint > cells( mesh_header.num_centities, 0 );

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all(
                     file, cells.data(), mesh_header.num_centities,
                     byte_offset + mesh_header.offsets[1] * sizeof( uint ),
                     mesh_header.disp[1] );
#else
    file.read( ( char * ) cells, ( mesh_header.num_centities ) * sizeof( uint ) );
#endif

    Array< uint > v;
    uint ci = 0;
    for ( uint i = 0; i < mesh_header.num_centities; i += mesh_header.num_entities )
    {
      v.clear();
      for ( uint j = 0; j < mesh_header.num_entities; ++j )
      {
        v.push_back( cells[i + j] );
      }
      editor.add_cell( ci++, &v[0] );
    }
  }

  if ( MPI::size() > 1 )
  {
    Array< uint > mapping( _mesh.size( 0 ) );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all(
                     file, mapping.data(), mesh_header.num_vertices,
                     byte_offset + mesh_header.offsets[2] * sizeof( uint ),
                     mesh_header.disp[2] );
#else
    file.read( ( char * ) mapping.data(), mesh_header.num_vertices * sizeof( uint ) );
#endif
    for ( VertexIterator v( _mesh ); !v.end(); ++v )
      _mesh.distdata()[0].set_map( v->index(), mapping[v->index()] );

    Array< uint > ghosts( 2 * mesh_header.num_ghosts );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all(
                     file, ghosts.data(), 2 * mesh_header.num_ghosts,
                     byte_offset + mesh_header.offsets[3] * sizeof( uint ),
                     mesh_header.disp[3] );
#else
    file.read( ( char * ) ghosts.data(), 2 * mesh_header.num_ghosts * sizeof( uint ) );
#endif
    for ( uint i = 0; i < 2 * mesh_header.num_ghosts; i += 2 )
    {
      _mesh.distdata()[0].set_ghost( ghosts[i], ghosts[i + 1] );
    }

    _mesh.distdata()[0].remap_shared_adj();
    _mesh.distdata()[0].finalize();
    _mesh.distdata()[0].valid_numbering = true;
  }

  editor.close();
  swap( mesh, _mesh );

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, LabelList< Function > & func )
{
  stream_t file = open_file( filename );

  // set the correct offset
  Checkpoint::offset_t byte_offset = chkp_header.offset_functions;

  FunctionsHeader functions_header;

  // read number of functions
#ifdef ENABLE_MPIIO
  byte_offset += MPI::file_read_at_all( file, &functions_header.count, 1,
                   byte_offset + MPI::rank() * sizeof( uint ), MPI::size() );
#else
  file.read( ( char * ) &functions_header.count, sizeof( uint ) );
#endif

  functions_header.offset.resize( functions_header.count );
  functions_header.names.resize( functions_header.count );

  // read all functions from the checkpoint
  for ( uint i = 0; i < functions_header.count; ++i )
  {
    functions_header.offset[i].resize( 3 );
    // read header and name
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, &functions_header.offset[i][0], 3,
                     byte_offset + MPI::rank() * 3 * sizeof( uint ),
                     MPI::size() * 3 );

    functions_header.names[i].resize( functions_header.name_length, '!' );
    byte_offset += MPI::file_read_at_all( file, &functions_header.names[i][0],
                     functions_header.name_length, ( byte_offset + MPI::rank()
                       * functions_header.name_length * sizeof( char ) ),
                     MPI::size() * functions_header.name_length );
#else
    file.read( ( char * ) &functions_header.offset[i][1], sizeof( uint ) );
#endif

    // find out if we 'requested' this function in the LabelList
    functions_header.names[i] = functions_header.names[i].substr(
                                  0, functions_header.names[i].find( '?' ) );
    uint index = func.size() + 1;
    for ( uint idx = 0; idx < func.size() ; ++idx )
    {
      if ( functions_header.names[i].size() == func[idx].second.size()
           and functions_header.names[i].find( func[idx].second ) == 0 )
      {
        index = idx;
        break;
      }
    }

    Array< real > values( functions_header.offset[i][1] );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, values.data(), values.size(),
                     byte_offset + functions_header.offset[i][0] * sizeof( real ),
                     functions_header.offset[i][2] );
#else
    file.read( ( char * ) values.data(), functions_header.offset[i][1] * sizeof( real ) );
#endif

    if ( index < func.size() )
    {
      if ( functions_header.offset[i][1] != func[index].first->vector().local_size() )
      {
        error( "Size mismatch while reloading functions from checkpoint:\n"
               "\tExpected : %d\n\tRead     : %d\n",
               func[index].first->vector().local_size(), functions_header.offset[i][1] );
      }

      func[index].first->vector().set( values.data() );
      func[index].first->vector().apply();
    }
  }

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, LabelList< GenericVector > & vec )
{
  stream_t file = open_file( filename );

  // set the correct offset
  Checkpoint::offset_t byte_offset = chkp_header.offset_vectors;

  VectorsHeader vectors_header;

  // read number of functions
#ifdef ENABLE_MPIIO
  byte_offset += MPI::file_read_at_all( file, &vectors_header.count, 1,
                   byte_offset + MPI::rank() * sizeof( uint ), MPI::size() );
#else
  file.read( ( char * ) &vectors_header.count, sizeof( uint ) );
#endif

  vectors_header.offset.resize( vectors_header.count );
  vectors_header.names.resize( vectors_header.count );

  // read all functions from the checkpoint
  for ( uint i = 0; i < vectors_header.count; ++i )
  {
    vectors_header.offset[i].resize( 3 );

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, vectors_header.offset[i].data(),
                     3, byte_offset + MPI::rank() * 3 * sizeof( uint ),
                     MPI::size() * 3 );

    vectors_header.names[i].resize( vectors_header.name_length, '!' );
    byte_offset += MPI::file_read_at_all( file, &vectors_header.names[i][0],
                     vectors_header.name_length, ( byte_offset + MPI::rank()
                       * vectors_header.name_length * sizeof( char ) ),
                     MPI::size() * vectors_header.name_length );
#else
    file.read( ( char * ) &vectors_header.offset[i][1], sizeof( uint ) );
#endif

    // find out if we 'requested' this function in the LabelList
    vectors_header.names[i] = vectors_header.names[i].substr(
                                0, vectors_header.names[i].find( '?' ) );
    uint index = vec.size() + 1;
    for ( uint idx = 0; idx < vec.size() ; ++idx )
    {
      if ( vectors_header.names[i].size() == vec[idx].second.size()
           and vectors_header.names[i].find( vec[idx].second ) == 0 )
      {
        index = idx;
        break;
      }
    }

    Array< real > values( vectors_header.offset[i][1] );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, values.data(), values.size(),
                     byte_offset + vectors_header.offset[i][0] * sizeof( real ),
                     vectors_header.offset[i][2] );
#else
    file.read( ( char * ) values.data(), vectors_header.offset[i][1] * sizeof( real ) );
#endif

    if ( index < vec.size() )
    {
      if ( vectors_header.offset[i][1] != vec[index].first->local_size() )
      {
        error( "Size mismatch while reloading vectors from checkpoint:\n"
               "\tExpected : %d\n\tRead     : %d\n",
               vec[index].first->local_size(), vectors_header.offset[i][1] );
      }

      vec[index].first->set( values.data() );
      vec[index].first->apply();
    }
  }

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::fill_headers( real const t,
                               Mesh & mesh, MeshHeader & mesh_header,
                               LabelList< Function > & func,
                               FunctionsHeader & functions_header,
                               LabelList< GenericVector > & vec,
                               VectorsHeader & vectors_header )
{
  // mesh header
  {
    mesh_header.num_coords    = mesh.size( 0 ) * mesh.geometry_dimension();
    mesh_header.num_entities  = mesh.type().num_entities( 0 );
    mesh_header.num_centities = mesh.num_cells() * mesh_header.num_entities;
    mesh_header.type          = mesh.type().cellType();
    mesh_header.tdim          = mesh.topology_dimension();
    mesh_header.gdim          = mesh.geometry_dimension();
    mesh_header.num_vertices  = mesh.size( 0 );
    mesh_header.num_cells     = mesh.num_cells();
    mesh_header.num_ghosts    = mesh.topology().num_ghost( 0 );

#ifdef ENABLE_MPIIO
    uint local_data[4] = { mesh_header.num_coords,
                           mesh_header.num_centities,
                           mesh_header.num_vertices,
                           ( 2 * mesh_header.num_ghosts ) };

    memset( &mesh_header.offsets[0], 0, 4 * sizeof( uint ) );
    MPI::exscan_sum( &local_data[0], &mesh_header.offsets[0], 4 );

    memset( &mesh_header.disp[0], 0, 4 * sizeof( uint ) );
    MPI::all_reduce< MPI::sum >( local_data, mesh_header.disp, 4 );
#endif
  }

  // functions header
  {
    functions_header.offset.resize( func.size() );
    functions_header.names.resize( func.size() );

    // fill header
    functions_header.count = func.size();
    functions_header.total_offset = MPI::size() * sizeof( uint );

    for ( uint i = 0; i < func.size(); ++i )
    {
      Function * f = func[i].first;
      functions_header.offset[i].resize( 3 );
      functions_header.offset[i][0] = f->vector().offset();
      functions_header.offset[i][1] = f->vector().local_size();
      functions_header.offset[i][2] = f->vector().size();
      functions_header.names.push_back( func[i].second );

      // header offset
      functions_header.total_offset += MPI::size() * ( 3 * sizeof( uint )
                              + functions_header.name_length * sizeof( char ) );

      // function offset
      functions_header.total_offset+= functions_header.offset[i][2]
                                       * sizeof( real );
    }
  }

  // vectors header
  {
    vectors_header.offset.resize( func.size() );
    vectors_header.names.resize( 0 );

    // fill header
    vectors_header.count = vec.size();
    vectors_header.total_offset = MPI::size() * sizeof( uint );

    for ( uint i = 0; i < vec.size(); ++i )
    {
      GenericVector * f = vec[i].first;
      vectors_header.offset[i].resize( 3 );
      vectors_header.offset[i][0] = f->offset();
      vectors_header.offset[i][1] = f->local_size();
      vectors_header.offset[i][2] = f->size();
      vectors_header.names.push_back( vec[i].second );

      // header offset
      vectors_header.total_offset += MPI::size() * ( 3 * sizeof( uint )
                                + vectors_header.name_length * sizeof( char ) );

      // vector offset
      vectors_header.total_offset += vectors_header.offset[i][2]
                                     * sizeof( real );
    }
  }

  // checkpoint header
  chkp_header.time = t;
  chkp_header.offset_psystem   = sizeof( Checkpoint::CheckpointHeader );
  chkp_header.offset_mesh     += chkp_header.offset_psystem + sizeof( uint );
  chkp_header.offset_functions = chkp_header.offset_mesh
                                 + MPI::size() * sizeof( MeshHeader )
                                 + mesh_header.disp[0] * sizeof( real )
                                 + mesh_header.disp[1] * sizeof( uint );
  if ( MPI::size() > 1 )
    chkp_header.offset_functions += mesh_header.disp[2] * sizeof( uint )
                                    + mesh_header.disp[3] * sizeof( uint );
  chkp_header.offset_vectors   = chkp_header.offset_functions
                                 + functions_header.total_offset;
}

//-----------------------------------------------------------------------------

void Checkpoint::write_mesh( stream_t file, offset_t & byte_offset,
                             Mesh & mesh, MeshHeader & mesh_header )
{
#ifdef ENABLE_MPIIO
  MPI::file_write_at_all( file, mesh_header,
                     byte_offset + MPI::rank() * sizeof( MeshHeader ) );
  byte_offset += MPI::size() * sizeof( MeshHeader );

  byte_offset += MPI::file_write_at_all(
                   file, mesh.geometry().coordinates(), mesh_header.num_coords,
                   byte_offset + mesh_header.offsets[0] * sizeof( real ),
                   mesh_header.disp[0] );

  Array< uint > cell_data;
  for ( uint c1 = 0; c1 < mesh.cells().size(); ++c1 )
    cell_data.append( mesh.cells()[c1].begin(), mesh.cells()[c1].end() );

  byte_offset += MPI::file_write_at_all(
                   file, cell_data.data(), mesh_header.num_centities,
                   byte_offset + mesh_header.offsets[1] * sizeof( uint ),
                   mesh_header.disp[1] );
#else
  file.write( static_cast< char * >( &mesh_header, sizeof( MeshHeader ) ) );
  file.write( static_cast< char * >( mesh.geometry().coordinates(),
             mesh_header.num_coords * sizeof( real ) ) );
  file.write( static_cast< char * >( mesh.cells(), mesh_header.num_centities * sizeof( uint ) ) );
#endif

  if ( MPI::size() > 1 )
  {
    Array< uint > mapping( mesh.size( 0 ) );
    for ( VertexIterator v( mesh ); !v.end(); ++v )
    {
      mapping[v->index()] = v->global_index();
    }

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all(
                     file, mapping.data(), mesh_header.num_vertices,
                     byte_offset + mesh_header.offsets[2] * sizeof( uint ),
                     mesh_header.disp[2] );
#else
    file.write( static_cast< char * >( mapping.data(), mesh_header.num_vertices * sizeof( uint ) ) );
#endif

    Array< uint > ghosts( 2 * mesh_header.num_ghosts );
    uint * gp     = &ghosts[0];
    for ( GhostIterator g( mesh.distdata()[0] ); g.valid(); ++g )
    {
      *gp++ = g.index();
      *gp++ = g.owner();
    }

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all(
                     file, ghosts.data(), 2 * mesh_header.num_ghosts,
                     byte_offset + mesh_header.offsets[3] * sizeof( uint ),
                     mesh_header.disp[3] );
#else
    file.write( static_cast< char * >( ghosts.data(), 2 * mesh_header.num_ghosts * sizeof( uint ) ) );
#endif
  }
}

//-----------------------------------------------------------------------------

void Checkpoint::write_func( stream_t file, offset_t & byte_offset,
                             LabelList< Function > & func,
                             FunctionsHeader & functions_header )
{
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all( file, &functions_header.count, 1,
                     byte_offset + MPI::rank() * sizeof( uint ), MPI::size() );
#else
    file.write( static_cast< char * >( &functions_header.count ), sizeof( uint ) );
#endif

  for ( uint i = 0; i < func.size(); ++i )
  {
    Array< real > values( functions_header.offset[i][1] );
    func[i].first->vector().get( values.data() );

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all( file, functions_header.offset[i].data(),
                     3, byte_offset + MPI::rank() * 3 * sizeof( uint ),
                     MPI::size() * 3 );

    func[i].second.resize( functions_header.name_length, '?' );
    byte_offset += MPI::file_write_at_all( file, func[i].second.c_str(), functions_header.name_length,
                     byte_offset + MPI::rank() * functions_header.name_length * sizeof( char ),
                     MPI::size() * functions_header.name_length );

    byte_offset += MPI::file_write_at_all( file, values.data(), values.size(),
                     byte_offset + functions_header.offset[i][0] * sizeof( real ),
                     func[i].first->vector().size() );
#else
    uint local_size = func[i].first->vector().local_size();
    file.write( static_cast< char * >( &local_size, sizeof( uint ) ) );
    file.write( static_cast< char * >( values.data(),
               func[i].first->vector().local_size() * sizeof( real ) ) );
#endif
  }
}

//-----------------------------------------------------------------------------

void Checkpoint::write_vect( stream_t file, offset_t & byte_offset,
                             LabelList< GenericVector > & vec,
                             VectorsHeader & vectors_header )
{
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all( file, &vectors_header.count, 1,
                     byte_offset + MPI::rank() * sizeof( uint ), MPI::size() );
#else
    file.write( static_cast< char * >( &vectors_header.count ), sizeof( uint ) );
#endif

  for ( uint i = 0; i < vec.size(); ++i )
  {
    Array< real > values( vectors_header.offset[i][1] );
    vec[i].first->get( values.data() );

#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_write_at_all( file, vectors_header.offset[i].data(),
                     3, byte_offset + MPI::rank() * 3 * sizeof( uint ),
                     MPI::size() * 3 );

    vec[i].second.resize( vectors_header.name_length, '?' );
    byte_offset += MPI::file_write_at_all(
                     file, vec[i].second.c_str(), vectors_header.name_length,
                     byte_offset + MPI::rank() * vectors_header.name_length * sizeof( char ),
                     MPI::size() * vectors_header.name_length );

    byte_offset += MPI::file_write_at_all( file, values.data(), values.size(),
                     byte_offset + vectors_header.offset[i][0] * sizeof( real ),
                     vectors_header.offset[i][2] );
#else
    uint local_size = it->first->local_size();
    file.write( static_cast< char * >( &local_size, sizeof( uint ) ) );
    file.write( static_cast< char * >( values.data(), vec[i].first->local_size() * sizeof( real ) ) );
#endif
  }
}

//-----------------------------------------------------------------------------

std::string Checkpoint::get_filename( std::string filebasename )
{
  // if the input filename contains the extension, remove it
  if ( filebasename.rfind( ".chkp" ) == filebasename.size() - 5 )
    filebasename.resize( filebasename.size() - 5 );

  // build filename
  std::ostringstream filename;
  filename << filebasename << n_;

#ifndef ENABLE_MPIIO
  if ( MPI::size() > 1 )
    filename << "_" << MPI::rank();
#endif

  filename << ".chkp";

  return filename.str();
}

//-----------------------------------------------------------------------------

Checkpoint::stream_t Checkpoint::open_file( std::string filebasename )
{
#ifdef ENABLE_MPIIO
  stream_t file = MPI_File();
  MPI::file_open( file, get_filename( filebasename ), MPI_MODE_RDONLY );

  // load header
  byte_offset = MPI::file_read_all(
                  file, chkp_header, sizeof( CheckpointHeader ) );
#else
  stream_t file = stream_t( get_filename( filebasename ), std::ofstream::binary );

  // load header
  file.read( static_cast< char * >( &chkp_header ),
               sizeof( CheckpointHeader ) );
#endif
  return file;
}

//-----------------------------------------------------------------------------

void Checkpoint::close_file( stream_t & file )
{
#ifdef ENABLE_MPIIO
  MPI::file_close( file );
#else
  file.close();
#endif
}

//-----------------------------------------------------------------------------

} // end namespace dolfin
