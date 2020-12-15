// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHECKPOINT_H
#define __DOLFIN_CHECKPOINT_H

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/CellType.h>

#include <fstream>
#include <string>

namespace dolfin
{

class Mesh;
class GenericVector;
class Function;

class Checkpoint
{
public:
  using MeshMap     = _ordered_map< std::string, Mesh * >;
  using FunctionMap = _ordered_map< std::string, Function * >;
  using VectorMap   = _ordered_map< std::string, GenericVector * >;

#ifdef ENABLE_MPIIO
  using stream_t = MPI_File;
  using offset_t = MPI_Offset;
#else
  using stream_t = std::ofstream;
  using offset_t = long long;
#endif

  static uint32_t const NAME_LENGTH = 256;

  struct CheckpointHeader
  {
    CheckpointHeader() = default;
    void disp() const;

    double   time{ 0.0 };
    uint32_t magic;
    uint32_t pe_size{ 0 };
    uint32_t num_meshes{ 0 };
    uint32_t num_functions{ 0 };
    uint32_t num_vectors{ 0 };
    offset_t offset_psystem{ 0 };
    offset_t offset_mesh{ 0 };
    offset_t offset_functions{ 0 };
    offset_t offset_vectors{ 0 };
  };

  struct MeshHeader
  {
    MeshHeader() = default;
    void disp() const;

    CellType::Type type{ CellType::point };
    uint32_t tdim{ 0 };
    uint32_t gdim{ 0 };
    uint32_t num_vertices{ 0 };
    uint32_t num_cells{ 0 };
    uint32_t num_entities{ 0 };
    uint32_t num_centities{ 0 };
    uint32_t num_coords{ 0 };
    uint32_t num_ghosts{ 0 };
    char     name[NAME_LENGTH];
  #ifdef ENABLE_MPIIO
    uint32_t offsets[4]{ 0, 0, 0, 0 };
    uint32_t displacement[4]{ 0, 0, 0, 0 };
  #endif
  };

  struct FunctionHeader
  {
    FunctionHeader() = default;
    void disp() const;

    uint32_t dim{ 0 };
    uint32_t size{ 0 };
    uint32_t offset[3]{ 0, 0, 0 };
    char     name[NAME_LENGTH];
  };

  struct VectorHeader
  {
    VectorHeader() = default;
    void disp() const;

    uint32_t offset[3]{ 0, 0, 0 };
    char     name[NAME_LENGTH];
  };

public:
  ///
  Checkpoint() = default;

  ///
  ~Checkpoint() = default;

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
  void increment_counter();

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
  uint n_{ 0 };

  CheckpointHeader        chkp_header;
  Array< MeshHeader >     mesh_header;
  Array< FunctionHeader > functions_header;
  Array< VectorHeader >   vectors_header;
};

}
#endif
