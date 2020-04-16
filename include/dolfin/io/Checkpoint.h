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
  typedef std::map< std::string, Function * >      FunctionMap;
  typedef std::map< std::string, GenericVector * > VectorMap;

#ifdef ENABLE_MPIIO
  typedef MPI_File   stream_t;
  typedef MPI_Offset offset_t;
#else
  typedef std::ofstream stream_t;
  typedef long long     offset_t;
#endif

  struct CheckpointHeader
  {
    real     time;
    offset_t offset_psystem;
    offset_t offset_mesh;
    offset_t offset_functions;
    offset_t offset_vectors;
  };

  struct MeshHeader
  {
    CellType::Type type;
    uint tdim;
    uint gdim;
    uint num_vertices;
    uint num_cells;
    uint num_entities;
    uint num_centities;
    uint num_coords;
    uint num_ghosts;
  #ifdef ENABLE_MPIIO
    uint offsets[4];
    uint disp[4];
  #endif
  };

  struct FunctionsHeader
  {
    uint                   count;
    Checkpoint::offset_t   total_offset;
    Array< Array< uint > > offset;
    Array< std::string >   names;

    static uint const name_length = 20;
  };

  typedef FunctionsHeader VectorsHeader;

public:
  ///
  Checkpoint();

  ///
  ~Checkpoint();

  ///
  void write( std::string filename, real const t, Mesh & mesh,
              FunctionMap & func, VectorMap & vec );

  ///
  void load_parametersystem( std::string filename );

  ///
  void load( std::string filename, Mesh & mesh );

  ///
  void load( std::string filename, FunctionMap & func );

  ///
  void load( std::string filename, VectorMap & vec );

  ///
  uint id() const;

  ///
  real restart_time() const;

private:
  void fill_headers( real const t, Mesh & mesh,
                     FunctionMap & func, VectorMap & vec );

  void write( stream_t file, offset_t & byte_offset, Mesh & mesh );

  void write( stream_t file, offset_t & byte_offset, FunctionMap & func );

  void write( stream_t file, offset_t & byte_offset, VectorMap & vec );

  std::string build_filename( std::string filename );
  stream_t    load_file( std::string & filename );
  void        close_file( stream_t & file );

private:
  uint n_;

  CheckpointHeader chkp_header;
  MeshHeader       mesh_header;
  FunctionsHeader  functions_header;
  VectorsHeader    vectors_header;
};

}
#endif
