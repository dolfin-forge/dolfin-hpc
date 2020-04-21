// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHECKPOINT_H
#define __DOLFIN_CHECKPOINT_H

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/CellType.h>

#include <fstream>
#include <map>
#include <string>

namespace dolfin
{

class Mesh;
class GenericVector;
class Function;

class Checkpoint
{
public:
  typedef std::map< std::string, Mesh * >          MeshMap;
  typedef std::map< std::string, Function * >      FunctionMap;
  typedef std::map< std::string, GenericVector * > VectorMap;

#ifdef ENABLE_MPIIO
  typedef MPI_File   stream_t;
  typedef MPI_Offset offset_t;
#else
  typedef std::ofstream stream_t;
  typedef long long     offset_t;
#endif

  static uint32_t const NAME_LENGTH = 20;

  struct CheckpointHeader
  {
    CheckpointHeader();
    void disp() const;

    double   time;
    uint32_t pe_size;
    uint32_t num_meshes;
    uint32_t num_functions;
    uint32_t num_vectors;
    offset_t offset_psystem;
    offset_t offset_mesh;
    offset_t offset_functions;
    offset_t offset_vectors;
  };

  struct MeshHeader
  {
    MeshHeader();
    void disp() const;

    CellType::Type type;
    uint32_t tdim;
    uint32_t gdim;
    uint32_t num_vertices;
    uint32_t num_cells;
    uint32_t num_entities;
    uint32_t num_centities;
    uint32_t num_coords;
    uint32_t num_ghosts;
    char     name[NAME_LENGTH];
  #ifdef ENABLE_MPIIO
    uint32_t offsets[4];
    uint32_t displacement[4];
  #endif
  };

  struct FunctionHeader
  {
    FunctionHeader();
    void disp() const;

    uint32_t dim;
    uint32_t size;
    uint32_t offset[3];
    char     name[NAME_LENGTH];
  };

  struct VectorHeader
  {
    VectorHeader();
    void disp() const;

    uint32_t offset[3];
    char     name[NAME_LENGTH];
  };

public:
  ///
  Checkpoint();

  ///
  ~Checkpoint();

  ///
  void write( std::string filename, real const t, MeshMap & meshes,
              FunctionMap & func, VectorMap & vec );

  ///
  void load_header( std::string filename );

  ///
  void load_parametersystem( std::string filename );

  ///
  void load( std::string filename, MeshMap const & meshes );

  ///
  void load( std::string filename, FunctionMap const & func );

  ///
  void load( std::string filename, VectorMap const & vec );

  ///
  real time() const;

  ///
  void reset_counter();

  ///
  CheckpointHeader        const & get_header() const;
  Array< MeshHeader >     const & get_mesh_header() const;
  Array< FunctionHeader > const & get_function_header() const;
  Array< VectorHeader >   const & get_vector_header() const;

private:
  void fill_headers( real const t, uint param_size, MeshMap & meshes,
                     FunctionMap & func, VectorMap & vec );

  void write( stream_t file, offset_t & byte_offset, MeshMap & meshes );

  void write( stream_t file, offset_t & byte_offset, FunctionMap & func );

  void write( stream_t file, offset_t & byte_offset, VectorMap & vec );

  std::string build_filename( std::string filename );
  stream_t    load_file( std::string & filename );
  void        close_file( stream_t & file );

private:
  uint n_;

  CheckpointHeader        chkp_header;
  Array< MeshHeader >     mesh_header;
  Array< FunctionHeader > functions_header;
  Array< VectorHeader >   vectors_header;
};

}
#endif
