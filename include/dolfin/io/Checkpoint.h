// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHECKPOINT_H
#define __DOLFIN_CHECKPOINT_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/celltypes/CellType.h>

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

    double   time { 0.0 };
    uint32_t magic;
    uint32_t pe_size { 0 };
    uint32_t num_meshes { 0 };
    uint32_t num_functions { 0 };
    uint32_t num_vectors { 0 };
    offset_t offset_psystem { 0 };
    offset_t offset_mesh { 0 };
    offset_t offset_functions { 0 };
    offset_t offset_vectors { 0 };
  };

  struct MeshHeader
  {
    MeshHeader() = default;
    void disp() const;

    CellType::Type type { CellType::point };
    uint32_t       tdim { 0 };
    uint32_t       gdim { 0 };
    uint32_t       num_vertices { 0 };
    uint32_t       num_cells { 0 };
    uint32_t       num_entities { 0 };
    uint32_t       num_centities { 0 };
    uint32_t       num_coords { 0 };
    uint32_t       num_ghosts { 0 };
    char           name[NAME_LENGTH];
#ifdef ENABLE_MPIIO
    uint32_t offsets[4] { 0, 0, 0, 0 };
    uint32_t displacement[4] { 0, 0, 0, 0 };
#endif
  };

  struct FunctionHeader
  {
    FunctionHeader() = default;
    void disp() const;

    uint32_t dim { 0 };
    uint32_t size { 0 };
    uint32_t offset[3] { 0, 0, 0 };
    char     name[NAME_LENGTH];
  };

  struct VectorHeader
  {
    VectorHeader() = default;
    void disp() const;

    uint32_t offset[3] { 0, 0, 0 };
    char     name[NAME_LENGTH];
  };

public:
  ///
  Checkpoint() = default;

  ///
  ~Checkpoint() = default;

  ///
  void write( std::string   filename,
              real const    t,
              MeshMap &     meshes,
              FunctionMap & func,
              VectorMap &   vec );

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
  auto time() const -> real;

  ///
  void reset_counter();

  ///
  void increment_counter();

  ///
  auto get_header() const -> CheckpointHeader const &;
  auto get_mesh_header() const -> std::vector< MeshHeader > const &;
  auto get_function_header() const -> std::vector< FunctionHeader > const &;
  auto get_vector_header() const -> std::vector< VectorHeader > const &;

private:
  void fill_headers( real const    t,
                     size_t        param_size,
                     MeshMap &     meshes,
                     FunctionMap & func,
                     VectorMap &   vec );

  void write( stream_t file, offset_t & byte_offset, MeshMap & meshes );
  void write( stream_t file, offset_t & byte_offset, FunctionMap & func );
  void write( stream_t file, offset_t & byte_offset, VectorMap & vec );

  auto build_filename( std::string filename ) -> std::string;
  auto load_file( std::string & filename ) -> stream_t;
  void close_file( stream_t & file );

private:
  size_t n_ { 0 };

  CheckpointHeader              chkp_header;
  std::vector< MeshHeader >     mesh_header;
  std::vector< FunctionHeader > functions_header;
  std::vector< VectorHeader >   vectors_header;
};

}
#endif
