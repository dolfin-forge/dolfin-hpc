// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHECKPOINT_H
#define __DOLFIN_CHECKPOINT_H

#include <dolfin/common/Array.h>
#include <dolfin/common/Label.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/CellType.h>

#include <fstream>
#include <string>

namespace dolfin
{

class Mesh;

class Checkpoint
{
public:
#ifdef ENABLE_MPIIO
  typedef MPI_File stream_t;
  typedef MPI_Offset offset_t;
#else
  typedef std::ofstream stream_t;
  typedef long long offset_t;
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
              LabelList< Function > & func, LabelList< GenericVector > & vec );

  ///
  void load_parametersystem( std::string filename );

  ///
  void load( std::string filename, Mesh & mesh );

  ///
  void load( std::string filename, LabelList< Function > & func );

  ///
  void load( std::string filename, LabelList< GenericVector > & vec );

  ///
  inline bool restart()
  {
    return true;//state_ == RESTART;
  }

  ///
  inline dolfin::real restart_time()
  {
    return chkp_header.time;
  }

private:
  void fill_headers( real const t,
                     Mesh & mesh, MeshHeader & mesh_header,
                     LabelList< Function > & func,
                     FunctionsHeader & functions_header,
                     LabelList< GenericVector > & vec,
                     VectorsHeader & vectors_header );

  void write_mesh( stream_t file, Checkpoint::offset_t & byte_offset,
                   Mesh & mesh,
                   MeshHeader & mesh_header );

  void write_func( stream_t file, Checkpoint::offset_t & byte_offset,
                   LabelList< Function > & func,
                   FunctionsHeader & functions_header );

  void write_vect( stream_t file, Checkpoint::offset_t & byte_offset,
                   LabelList< GenericVector > & vec,
                   VectorsHeader & vectors_header );

  std::string get_filename( std::string filebasename );
  stream_t open_file( std::string filebasename );
  void close_file( stream_t & file );

private:
  uint n_;

  offset_t byte_offset;
  CheckpointHeader chkp_header;
};

}
#endif
