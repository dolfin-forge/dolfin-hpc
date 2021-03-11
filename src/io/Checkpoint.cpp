// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/Checkpoint.h>

#include <dolfin/common/timing.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>
#include <dolfin/parameter/parameters.h>

#include <sstream>

uint32_t const CHECKPOINT_MAGIC = 0xDEADC0DE;

#ifdef DOLFIN_HAVE_MPI

namespace dolfin
{

/* Checkpoint file format:
 * 1. CheckpointHeader
 * 2. ParameterSystem
 * 2a. ParameterSystem size (uint)
 * 2b. ParameterSystem (string)
 * 3. Meshes (for each mesh)
 * 3a. MeshHeader
 * 3b. Mesh
 * 4. Functions (for each Function)
 * 4a. FunctionHeader
 * 4b. Function
 * 5. Vectors (for each Vector)
 * 5a. VectorHeader
 * 5b. Vector
 */

//-----------------------------------------------------------------------------

void Checkpoint::write( std::string filename, real const t, MeshMap & meshes,
                        FunctionMap & func, VectorMap & vec )
{
  tic();

  dolfin_set( "checkpoint_id", n_ );
  dolfin_set( "checkpoint_time", t );

  filename = build_filename( filename );
  message( "Writing checkpoint (name = %s, id = %d) at time %g",
           filename.c_str(), n_ - 1, t );

  // deserialize ParameterSystem
  std::string parameters  = ParameterSystem::parameters.serialize();
  uint        param_size  = parameters.size() * sizeof( char );

  fill_headers( t, param_size, meshes, func, vec );

#ifndef ENABLE_MPIIO
  error("Checkpointing is only implemented with MPI I/O enabled.");
  // stream_t file = stream_t( filename, std::ofstream::binary );

  // // write header
  // file.write( static_cast< char * >( &chkp_header ), sizeof( CheckpointHeader ) );

  // // write ParameterSystem
  // file.write( static_cast< char * >( &param_size ), sizeof( uint ) );
  // file.write( static_cast< char * >( parameters.c_str() ), param_size );
#else
  stream_t file = MPI_File();
  MPI::file_open( file, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE );

  // write header
  offset_t byte_offset = MPI::file_write_all( file, chkp_header,
                                              sizeof( CheckpointHeader ) );

  dolfin_assert( byte_offset == chkp_header.offset_psystem );

  // write ParameterSystem
  byte_offset += MPI::file_write_all( file, param_size );
  byte_offset += MPI::file_write_all( file, parameters.c_str()[0], param_size );
#endif

  write( file, byte_offset, meshes );
  write( file, byte_offset, func );
  write( file, byte_offset, vec );

  close_file( file );

  message( "Writing checkpoint (%s) took %f seconds.", filename.c_str(), toc() );
}

//-----------------------------------------------------------------------------

void Checkpoint::load_header( std::string filename )
{
  stream_t file = load_file( filename );

  // load mesh headers
  {
    offset_t byte_offset = chkp_header.offset_mesh;
    mesh_header.resize( chkp_header.num_meshes );

    for ( MeshHeader & mesh_hdr : mesh_header )
    {
      MPI::file_read_at_all( file, mesh_hdr, byte_offset );
      byte_offset += chkp_header.pe_size * sizeof( MeshHeader );
      byte_offset += sizeof( real ) * mesh_hdr.displacement[0];
      byte_offset += sizeof( uint ) * mesh_hdr.displacement[1];
      if ( chkp_header.pe_size > 1 )
      {
        byte_offset += sizeof( uint ) * mesh_hdr.displacement[2];
        byte_offset += sizeof( uint ) * mesh_hdr.displacement[3];
      }
    }
  }

  // load function headers
  {
    offset_t byte_offset  = chkp_header.offset_functions;
    functions_header.resize( chkp_header.num_functions );

    for ( FunctionHeader & function_hdr : functions_header )
    {
      MPI::file_read_at_all( file, function_hdr, byte_offset );
      byte_offset += chkp_header.pe_size * sizeof( FunctionHeader );
      byte_offset += sizeof( real ) * function_hdr.offset[2];
    }
  }

  // load vector headers
  {
    offset_t byte_offset  = chkp_header.offset_vectors;
    vectors_header.resize( chkp_header.num_vectors );

    for ( VectorHeader & vector_hdr : vectors_header )
    {
      MPI::file_read_at_all( file, vector_hdr, byte_offset );
      byte_offset += MPI::size() * sizeof( VectorHeader );
      byte_offset += sizeof( real ) * vector_hdr.offset[2];
    }
  }

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load_parametersystem( std::string filename )
{
  stream_t file = load_file( filename );

  offset_t byte_offset = chkp_header.offset_psystem;
  uint     param_size  = 0;

#ifdef ENABLE_MPIIO
  // load ParameterSystem
  byte_offset += MPI::file_read_at_all( file, param_size, byte_offset );
  std::vector< char > p( param_size / sizeof( char ) );
  byte_offset += MPI::file_read_at_all( file, &p[0],
                                        param_size, byte_offset, param_size );
#else
  error("Checkpointing is only implemented with MPI I/O enabled.");
  // load ParameterSystem
  // file.read( static_cast< char * >( &param_size ), sizeof( uint ) );
  // std::vector< char > p( param_size / sizeof( char ) );
  // file.read( static_cast< char * >( p.data() ), param_size );
#endif

  ParameterSystem::parameters.deserialize( std::string( p.data(), p.size() ) );
  message( 1, "Checkpoint: Loaded ParameterSystem from time %g", chkp_header.time );

  n_ = dolfin_get<size_t>( "checkpoint_id" );

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, MeshMap const & meshes )
{
  stream_t file = load_file( filename );

  offset_t byte_offset = chkp_header.offset_mesh;
  size_t     loaded_count = 0;
  mesh_header.resize( chkp_header.num_meshes );

  if ( chkp_header.pe_size != MPI::size() )
    error( "Checkpoint::load(Mesh): PE size does not match (%u != %u)",
           chkp_header.pe_size, MPI::size() );

  for ( MeshHeader & mesh_hdr : mesh_header )
  {
    // load header
#ifdef ENABLE_MPIIO
    MPI::file_read_at_all( file, mesh_hdr,
                           byte_offset + MPI::rank() * sizeof( MeshHeader ) );
    byte_offset += MPI::size() * sizeof( MeshHeader );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // file.read( ( char * ) &mesh_hdr, sizeof( MeshHeader ) );
#endif

    // find out if we 'requested' this function in the FunctionMap
    std::string name( mesh_hdr.name, NAME_LENGTH );
    name = name.substr( 0, name.find( '?' ) );
    MeshMap::const_iterator m_ = meshes.find( name );

    if ( m_ != meshes.end() )
    {
      Mesh _mesh;
      MeshEditor editor( _mesh, mesh_hdr.type, mesh_hdr.gdim );
      editor.init_vertices( mesh_hdr.num_vertices );

      // load coords
      {
        std::vector< real > coords( mesh_hdr.num_coords, 0.0 );

#ifdef ENABLE_MPIIO
        byte_offset += MPI::file_read_at_all(
                         file, coords.data(), mesh_hdr.num_coords,
                         byte_offset + mesh_hdr.offsets[0] * sizeof( real ),
                         mesh_hdr.displacement[0] );
 #else
        // file.read( ( char * ) coords, ( mesh_hdr.num_coords ) * sizeof( real ) );
 #endif

        for ( size_t i = 0; i < mesh_hdr.num_coords; i += mesh_hdr.gdim )
        {
          editor.add_vertex( i / mesh_hdr.gdim, coords.data() + i );
        }
      }

      // load cells
      {
        editor.init_cells( mesh_hdr.num_cells );
        std::vector< uint > cells( mesh_hdr.num_centities, 0 );

#ifdef ENABLE_MPIIO
        byte_offset += MPI::file_read_at_all(
                         file, cells.data(), mesh_hdr.num_centities,
                         byte_offset + mesh_hdr.offsets[1] * sizeof( uint ),
                         mesh_hdr.displacement[1] );
#else
       // file.read( ( char * ) cells, ( mesh_hdr.nuchkp_headerchkp_headerm_centities ) * sizeof( uint ) );
#endif

        std::vector< size_t > v;
        size_t ci = 0;
        for ( size_t i = 0; i < mesh_hdr.num_centities; i += mesh_hdr.num_entities )
        {
          v.clear();
          for ( size_t j = 0; j < mesh_hdr.num_entities; ++j )
          {
            v.push_back( cells[i + j] );
          }
          editor.add_cell( ci++, &v[0] );
        }
      }

      if ( MPI::size() > 1 )
      {
        std::vector< uint > mapping( _mesh.size( 0 ) );
#ifdef ENABLE_MPIIO
        byte_offset += MPI::file_read_at_all(
                         file, mapping.data(), mesh_hdr.num_vertices,
                         byte_offset + mesh_hdr.offsets[2] * sizeof( uint ),
                         mesh_hdr.displacement[2] );
#else
       // file.read( ( char * ) mapping.data(), mesh_hdr.num_vertices * sizeof( uint ) );
#endif
        for ( VertexIterator v( _mesh ); !v.end(); ++v )
          _mesh.distdata()[0].set_map( v->index(), mapping[v->index()] );

        std::vector< uint > ghosts( 2 * mesh_hdr.num_ghosts );
#ifdef ENABLE_MPIIO
        byte_offset += MPI::file_read_at_all(
                         file, ghosts.data(), 2 * mesh_hdr.num_ghosts,
                         byte_offset + mesh_hdr.offsets[3] * sizeof( uint ),
                         mesh_hdr.displacement[3] );
#else
       // file.read( ( char * ) ghosts.data(), 2 * mesh_hdr.num_ghosts * sizeof( uint ) );
#endif
        for ( size_t i = 0; i < 2 * mesh_hdr.num_ghosts; i += 2 )
        {
          _mesh.distdata()[0].set_ghost( ghosts[i], ghosts[i + 1] );
        }

        _mesh.distdata()[0].remap_shared_adj();
        _mesh.distdata()[0].finalize();
        _mesh.distdata()[0].valid_numbering = true;
      }

      editor.close();
      swap( *( m_->second ), _mesh );
      ++loaded_count;
    }
    else
    {
      // skip this mesh
      byte_offset += sizeof( real ) * mesh_hdr.displacement[0];
      byte_offset += sizeof( uint ) * mesh_hdr.displacement[1];
      if ( MPI::size() > 1 )
      {
        byte_offset += sizeof( uint ) * mesh_hdr.displacement[2];
        byte_offset += sizeof( uint ) * mesh_hdr.displacement[3];
      }
    }
  }

  message( 1, "Checkpoint: Loaded %d Mesh(es) from time %g",
              loaded_count, chkp_header.time );

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, FunctionMap const & func )
{
  stream_t file = load_file( filename );

  if ( chkp_header.pe_size != MPI::size() )
    error( "Checkpoint::load(Function): PE size does not match (%u != %u)",
           chkp_header.pe_size, MPI::size() );

  // set the correct offset
  offset_t byte_offset  = chkp_header.offset_functions;
  uint     loaded_count = 0;
  functions_header.resize( chkp_header.num_functions );

  // read all functions from the checkpoint
  for ( FunctionHeader & function_hdr : functions_header )
  {
#ifdef ENABLE_MPIIO
    // load header
    MPI::file_read_at_all( file, function_hdr,
                           byte_offset + MPI::rank() * sizeof( FunctionHeader ) );
    byte_offset += MPI::size() * sizeof( FunctionHeader );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // file.read( ( char * ) &function_hdr.offset[1], sizeof( uint ) );
#endif

    // load data
    std::vector< real > values( function_hdr.offset[1] );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, values.data(), values.size(),
                     byte_offset + function_hdr.offset[0] * sizeof( real ),
                     function_hdr.offset[2] );
#else
    // file.read( ( char * ) values.data(), function_hdr.offset[1] * sizeof( real ) );
#endif

    // find out if we 'requested' this function in the FunctionMap
    std::string name( function_hdr.name, NAME_LENGTH );
    name = name.substr( 0, name.find( '?' ) );
    FunctionMap::const_iterator f = func.find( name );

    if ( f != func.end() )
    {
      dolfin_assert( function_hdr.dim == f->second->value_size() );

      f->second->vector().set( values.data() );
      f->second->vector().apply();

      ++loaded_count;
    }
  }

  message( 1, "Checkpoint: Loaded %d Function(s) from time %g",
              loaded_count, chkp_header.time );

  close_file( file );
}

//-----------------------------------------------------------------------------

void Checkpoint::load( std::string filename, VectorMap const & vec )
{
  stream_t file = load_file( filename );

  if ( chkp_header.pe_size != MPI::size() )
    error( "Checkpoint::load(Vector): PE size does not match (%u != %u)",
           chkp_header.pe_size, MPI::size() );

  // set the correct offset
  offset_t byte_offset  = chkp_header.offset_vectors;
  uint     loaded_count = 0;
  vectors_header.resize( chkp_header.num_vectors );

  // read all functions from the checkpoint
  for ( VectorHeader & vector_hdr : vectors_header )
  {
#ifdef ENABLE_MPIIO
    // load header
    MPI::file_read_at_all( file, vector_hdr,
                           byte_offset + MPI::rank() * sizeof( VectorHeader ) );
    byte_offset += MPI::size() * sizeof( VectorHeader );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // file.read( ( char * ) &vector_hdrfset[i][1], sizeof( uint ) );
#endif

    // load data
    std::vector< real > values( vector_hdr.offset[1] );
#ifdef ENABLE_MPIIO
    byte_offset += MPI::file_read_at_all( file, values.data(), values.size(),
                     byte_offset + vector_hdr.offset[0] * sizeof( real ),
                     vector_hdr.offset[2] );
#else
    // file.read( ( char * ) values.data(), vector_hdr.offset[1] * sizeof( real ) );
#endif

    // find out if we 'requested' this function in the VectorMap
    std::string name( vector_hdr.name, NAME_LENGTH );
    name = name.substr( 0, name.find( '?' ) );
    VectorMap::const_iterator v = vec.find( name );

    if ( v != vec.end() )
    {
      if ( vector_hdr.offset[1] != v->second->local_size() )
      {
        error( "Size mismatch while reloading vectors from checkpoint:\n"
               "\tExpected : %d\n\tRead     : %d\n",
               v->second->local_size(), vector_hdr.offset[1] );
      }

      v->second->set( values.data() );
      v->second->apply();

      ++loaded_count;
    }
  }

  message( 1, "Checkpoint: Loaded %d Vector(s) from time %g",
              loaded_count, chkp_header.time );

  close_file( file );
}

//-----------------------------------------------------------------------------

auto Checkpoint::time() const -> real
{
  return chkp_header.time;
}

//-----------------------------------------------------------------------------

void Checkpoint::reset_counter()
{
  n_ = 0;
}

//-----------------------------------------------------------------------------

void Checkpoint::increment_counter()
{
  ++n_;
  message( 1, "Checkpoint: incremented checkpoint id to: %d", n_ );
}

//-----------------------------------------------------------------------------

auto Checkpoint::get_header() const -> Checkpoint::CheckpointHeader const &
{
  return chkp_header;
}

//-----------------------------------------------------------------------------

auto Checkpoint::get_mesh_header() const -> std::vector< Checkpoint::MeshHeader > const &
{
  return mesh_header;
}

//-----------------------------------------------------------------------------

auto Checkpoint::get_function_header() const -> std::vector< Checkpoint::FunctionHeader > const &
{
  return functions_header;
}

//-----------------------------------------------------------------------------

auto Checkpoint::get_vector_header() const -> std::vector< Checkpoint::VectorHeader > const &
{
  return vectors_header;
}

//-----------------------------------------------------------------------------

void Checkpoint::fill_headers( real const t, size_t param_size, MeshMap & meshes,
                               FunctionMap & func, VectorMap & vec )
{
  // mesh header
  mesh_header.resize( 0 );
  chkp_header.offset_mesh = 0;

  for ( MeshMap::iterator m = meshes.begin(); m != meshes.end(); ++m )
  {
    Mesh & mesh = *( m->second );
    MeshHeader hdr;

    hdr.num_coords    = mesh.size( 0 ) * mesh.geometry_dimension();
    hdr.num_entities  = mesh.type().num_entities( 0 );
    hdr.num_centities = mesh.num_cells() * hdr.num_entities;
    hdr.type          = mesh.type().cellType();
    hdr.tdim          = mesh.topology_dimension();
    hdr.gdim          = mesh.geometry_dimension();
    hdr.num_vertices  = mesh.size( 0 );
    hdr.num_cells     = mesh.num_cells();
    hdr.num_ghosts    = mesh.topology().num_ghost( 0 );

    std::string name( m->first );
    name.resize( NAME_LENGTH, '?' );
    name.copy( hdr.name, sizeof( hdr.name ) );

#ifdef ENABLE_MPIIO
    uint32_t local_data[4] = { hdr.num_coords, hdr.num_centities,
                               hdr.num_vertices, 2 * hdr.num_ghosts };

    MPI::exscan_sum( &local_data[0], &hdr.offsets[0], 4 );
    MPI::all_reduce< MPI::sum >( local_data, hdr.displacement, 4 );
#endif

    mesh_header.push_back( hdr );

    // mesh offset
    chkp_header.offset_mesh += MPI::size() * sizeof( MeshHeader )
                                + hdr.displacement[0] * sizeof( real )
                                + hdr.displacement[1] * sizeof( uint );

    if ( MPI::size() > 1 )
      chkp_header.offset_mesh += hdr.displacement[2] * sizeof( uint )
                                  + hdr.displacement[3] * sizeof( uint );
  }

  // fill functions header
  functions_header.resize( 0 );
  chkp_header.offset_functions = 0;

  for ( FunctionMap::iterator f = func.begin(); f != func.end(); ++f )
  {
    Function * fun = f->second;
    FunctionHeader hdr;

    hdr.dim       = fun->value_size();
    hdr.size      = fun->value_size() * fun->mesh().global_size(0);
    hdr.offset[0] = fun->vector().offset();
    hdr.offset[1] = fun->vector().local_size();
    hdr.offset[2] = fun->vector().size();

    std::string name( f->first );
    name.resize( NAME_LENGTH, '?' );
    name.copy( hdr.name, sizeof( hdr.name ) );

    functions_header.push_back( hdr );

    // function offset
    chkp_header.offset_functions += MPI::size() * sizeof( FunctionHeader );
    chkp_header.offset_functions += hdr.offset[2] * sizeof( real );
  }

  // fill functions header
  vectors_header.resize( 0 );
  chkp_header.offset_vectors = 0;

  for ( VectorMap::iterator v = vec.begin(); v != vec.end(); ++v )
  {
    GenericVector * vec = v->second;
    VectorHeader hdr;

    hdr.offset[0] = vec->offset();
    hdr.offset[1] = vec->local_size();
    hdr.offset[2] = vec->size();

    std::string name( v->first );
    name.resize( NAME_LENGTH, '?' );
    name.copy( hdr.name, sizeof( hdr.name ) );

    vectors_header.push_back( hdr );
  }

  // checkpoint header
  chkp_header.time          = t;
  chkp_header.magic         = CHECKPOINT_MAGIC;
  chkp_header.pe_size       = MPI::size();
  chkp_header.num_meshes    = meshes.size();
  chkp_header.num_functions = func.size();
  chkp_header.num_vectors   = vec.size();

  dolfin_assert( chkp_header.num_meshes    == mesh_header.size() );
  dolfin_assert( chkp_header.num_functions == functions_header.size() );
  dolfin_assert( chkp_header.num_vectors   == vectors_header.size() );

  chkp_header.offset_psystem   = sizeof( CheckpointHeader );
  chkp_header.offset_vectors   = chkp_header.offset_functions;
  chkp_header.offset_functions = chkp_header.offset_mesh;
  chkp_header.offset_mesh      = chkp_header.offset_psystem;
  chkp_header.offset_mesh      += param_size + sizeof( uint );
  chkp_header.offset_functions += chkp_header.offset_mesh;
  chkp_header.offset_vectors   += chkp_header.offset_functions;

  dolfin_assert( chkp_header.offset_psystem   < chkp_header.offset_mesh );
  dolfin_assert( chkp_header.offset_mesh      <= chkp_header.offset_functions );
  dolfin_assert( chkp_header.offset_functions <= chkp_header.offset_vectors );
}

//-----------------------------------------------------------------------------

void Checkpoint::write( stream_t file, offset_t & byte_offset, MeshMap & meshes )
{
  dolfin_assert( byte_offset == chkp_header.offset_mesh );

#ifdef ENABLE_MPIIO
  uint i = 0;
  for ( MeshMap::iterator m = meshes.begin(); m != meshes.end(); ++m, ++i )
  {
    MeshHeader & hdr = mesh_header[i];
    Mesh & mesh = *( m->second );

    MPI::file_write_at_all( file, hdr,
                            byte_offset + MPI::rank() * sizeof( MeshHeader ) );
    byte_offset += MPI::size() * sizeof( MeshHeader );

    byte_offset += MPI::file_write_at_all( file, mesh.geometry().coordinates(),
                     hdr.num_coords, byte_offset + hdr.offsets[0] * sizeof( real ),
                     hdr.displacement[0] );

    std::vector< uint > cell_data;
    for ( uint c1 = 0; c1 < mesh.cells().size(); ++c1 )
      append( cell_data, mesh.cells()[c1].begin(), mesh.cells()[c1].end() );

    byte_offset += MPI::file_write_at_all( file, cell_data.data(),
                     hdr.num_centities,
                     byte_offset + hdr.offsets[1] * sizeof( uint ), hdr.displacement[1] );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // file.write( static_cast< char * >( &hdr, sizeof( MeshHeader ) ) );
    // file.write( static_cast< char * >( mesh.geometry().coordinates(),
    //            hdr.num_coords * sizeof( real ) ) );
    // file.write( static_cast< char * >( mesh.cells(), hdr.num_centities * sizeof( uint ) ) );
#endif

    if ( MPI::size() > 1 )
    {
      std::vector< uint > mapping( mesh.size( 0 ) );
      for ( VertexIterator v( mesh ); !v.end(); ++v )
      {
        mapping[v->index()] = v->global_index();
      }

#ifdef ENABLE_MPIIO
      byte_offset += MPI::file_write_at_all( file, mapping.data(), hdr.num_vertices,
                       byte_offset + hdr.offsets[2] * sizeof( uint ), hdr.displacement[2] );
#else
      // file.write( static_cast< char * >( mapping.data(), hdr.num_vertices * sizeof( uint ) ) );
#endif

      std::vector< uint > ghosts( 2 * hdr.num_ghosts );
      uint * gp     = &ghosts[0];
      for ( GhostIterator g( mesh.distdata()[0] ); g.valid(); ++g )
      {
        *gp++ = g.index();
        *gp++ = g.owner();
      }

#ifdef ENABLE_MPIIO
      byte_offset += MPI::file_write_at_all( file, ghosts.data(), 2 * hdr.num_ghosts,
                       byte_offset + hdr.offsets[3] * sizeof( uint ), hdr.displacement[3] );
#else
      // file.write( static_cast< char * >( ghosts.data(), 2 * hdr.num_ghosts * sizeof( uint ) ) );
#endif
    }
  }
}

//-----------------------------------------------------------------------------

void Checkpoint::write( stream_t file, offset_t & byte_offset,
                        FunctionMap & func )
{
  dolfin_assert( byte_offset == chkp_header.offset_functions );

  uint i = 0;
  for ( FunctionMap::iterator f = func.begin(); f != func.end(); ++f, ++i )
  {
    dolfin_assert( i < chkp_header.num_functions );

    std::vector< real > values( functions_header[i].offset[1] );
    f->second->vector().get( values.data() );

#ifdef ENABLE_MPIIO
    // write header
    MPI::file_write_at_all( file, functions_header[i],
                            byte_offset + MPI::rank() * sizeof( FunctionHeader ) );
    byte_offset += MPI::size() * sizeof( FunctionHeader );

    // write data
    byte_offset += MPI::file_write_at_all( file, values.data(), values.size(),
                     byte_offset + functions_header[i].offset[0] * sizeof( real ),
                     functions_header[i].offset[2] );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // uint local_size = f->second->vector().local_size();
    // file.write( static_cast< char * >( &local_size, sizeof( uint ) ) );
    // file.write( static_cast< char * >( values.data(),
    //            f->second->vector().local_size() * sizeof( real ) ) );
#endif
  }
}

//-----------------------------------------------------------------------------

void Checkpoint::write( stream_t file, offset_t & byte_offset, VectorMap & vec )
{
  dolfin_assert( byte_offset == chkp_header.offset_vectors );

  uint i = 0;
  for ( VectorMap::iterator v = vec.begin(); v != vec.end(); ++v, ++i )
  {
    dolfin_assert( i < chkp_header.num_vectors );

    std::vector< real > values( vectors_header[i].offset[1] );
    v->second->get( values.data() );

#ifdef ENABLE_MPIIO
    // write header
    MPI::file_write_at_all( file, vectors_header[i],
                            byte_offset + MPI::rank() * sizeof( VectorHeader ) );
    byte_offset += MPI::size() * sizeof( VectorHeader );

    // write data
    byte_offset += MPI::file_write_at_all( file, values.data(), values.size(),
                     byte_offset + vectors_header[i].offset[0] * sizeof( real ),
                     vectors_header[i].offset[2] );
#else
    error("Checkpointing is only implemented with MPI I/O enabled.");
    // uint local_size = it->first->local_size();
    // file.write( static_cast< char * >( &local_size, sizeof( uint ) ) );
    // file.write( static_cast< char * >( values.data(), v->second->local_size() * sizeof( real ) ) );
#endif
  }
}

//-----------------------------------------------------------------------------

auto Checkpoint::build_filename( std::string filename ) -> std::string
{
  // if the input filename contains the extension, remove it
  if ( filename.size() > 5
       and filename.rfind( ".chkp" ) == filename.size() - 5 )
    filename.resize( filename.size() - 5 );

  // build filename
  std::ostringstream filename_;
  filename_ << filename << n_++;

#ifndef ENABLE_MPIIO
  error("Checkpointing is only implemented with MPI I/O enabled.");
  // if ( MPI::size() > 1 )
  //   filename_ << "_" << MPI::rank();
#endif

  filename_ << ".chkp";

  return filename_.str();
}

//-----------------------------------------------------------------------------

auto Checkpoint::load_file( std::string & filename ) -> Checkpoint::stream_t
{
  // if the input filename contains the extension, remove it
  if ( filename.size() > 5
       and filename.rfind( ".chkp" ) != filename.size() - 5 )
    filename = filename + ".chkp";

#ifdef ENABLE_MPIIO
  stream_t file = MPI_File();
  MPI::file_open( file, filename, MPI_MODE_RDONLY );

  // load header
  MPI::file_read_all( file, chkp_header, sizeof( CheckpointHeader ) );

  if ( chkp_header.magic != CHECKPOINT_MAGIC )
    error( "File \"%s\" does not seem to contain a checkpoint!", filename.c_str() );
#else
  error("Checkpointing is only implemented with MPI I/O enabled.");
  // stream_t file = stream_t( filename, std::ifstream::binary );

  // // load header
  // file.read( static_cast< char * >( &chkp_header ), sizeof( CheckpointHeader ) );
#endif

  return file;
}

//-----------------------------------------------------------------------------

void Checkpoint::close_file( stream_t & file )
{
#ifdef ENABLE_MPIIO
  MPI::file_close( file );
#else
  error("Checkpointing is only implemented with MPI I/O enabled.");
  // file.close();
#endif
}

//-----------------------------------------------------------------------------

void Checkpoint::CheckpointHeader::disp() const
{
  begin("");
  header( "CheckpointHeader" );
  message( "time:             %e", time );
  message( "magic:            %X", magic );
  message( "pe_size:          %u", pe_size );
  message( "num_meshes:       %u", num_meshes );
  message( "num_functions:    %u", num_functions );
  message( "num_vectors:      %u", num_vectors );
  message( "offset_psystem:   %d", offset_psystem );
  message( "offset_mesh:      %d", offset_mesh );
  message( "offset_functions: %d", offset_functions );
  message( "offset_vectors:   %d", offset_vectors );
  end();
}

//-----------------------------------------------------------------------------

void Checkpoint::MeshHeader::disp() const
{
  std::string name_( name, NAME_LENGTH );
  begin("");
  header( "MeshHeader" );
  message( "type:          %u", type );
  message( "tdim:          %u", tdim );
  message( "gdim:          %u", gdim );
  message( "num_vertices:  %u", num_vertices );
  message( "num_cells:     %u", num_cells );
  message( "num_entities:  %u", num_entities );
  message( "num_centities: %u", num_centities );
  message( "num_coords:    %u", num_coords );
  message( "num_ghosts:    %u", num_ghosts );
  message( "name:          %s", name_.substr( 0, name_.find( '?' ) ).c_str() );
#ifdef ENABLE_MPIIO
  message( "offset: {%u,%u,%u,%u}",
           offsets[0], offsets[1], offsets[2], offsets[3] );
  message( "disp:   {%u,%u,%u,%u}",
           displacement[0], displacement[1], displacement[2], displacement[3] );
#endif
  end();
}

//-----------------------------------------------------------------------------

void Checkpoint::FunctionHeader::disp() const
{
  std::string name_( name, NAME_LENGTH );
  begin("");
  header( "FunctionHeader" );
  message( "dim:    %u", dim );
  message( "size:   %u", size );
  message( "offset: {%u,%u,%u}", offset[0], offset[1], offset[2] );
  message( "name:   %s", name_.substr( 0, name_.find( '?' ) ).c_str() );
  end();
}

//-----------------------------------------------------------------------------

void Checkpoint::VectorHeader::disp() const
{
  std::string name_( name, NAME_LENGTH );
  begin("");
  header( "VectorHeader" );
  message( "offset: {%u,%u,%u}", offset[0], offset[1], offset[2] );
  message( "name:   %s", name_.substr( 0, name_.find( '?' ) ).c_str() );
  end();
}

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif // DOLFIN_HAVE_MPI
