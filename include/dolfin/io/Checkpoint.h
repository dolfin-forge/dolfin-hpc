// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-08
// Last changed: 2011-04-19

#ifndef __CHECKPOINT_H
#define __CHECKPOINT_H

#include <fstream>
#include <string>
#include <vector>

#include <dolfin/la/Vector.h>
#include <dolfin/function/Function.h>

#ifdef ENABLE_MPIIO
#include <mpi.h>
#endif

namespace dolfin
{

class Mesh;

class Checkpoint
{

public:

  ///
  Checkpoint();

  ///
  ~Checkpoint();

  ///
  void write(std::string fname, uint id, real t, Mesh& mesh,
             std::vector<Function *> func, std::vector<Vector *> vec,
             bool static_mesh = false);

  ///
  void restart(std::string fname);

  ///
  void load(Mesh& mesh);

  ///
  void load(std::vector<Function *> func);

  ///
  void load(std::vector<Vector *> vec);

  ///
  inline bool restart()
  {
    return state_ == RESTART;
  }

  ///
  inline dolfin::uint id()
  {
    if (state_ != RESTART)
    {
      error("Shut her down, Scotty, she's sucking mud again!");
    }
    return id_;
  }

  ///
  inline dolfin::real restart_time()
  {
    if (state_ != RESTART)
    {
      error("Shut her down, Scotty, she's sucking mud again!");
    }
    return t_;
  }

  ///
  inline void reset()
  {
    state_ = CHECKPOINT;
    restart_state_ = OPEN;
    hdr_initialized_ = false;
    disp_initialized_ = false;
  }

private:

#ifdef ENABLE_MPIIO
  typedef MPI_File chkp_outstream;
#else
  typedef std::ofstream chkp_outstream;
#endif

  void hdr_init(Mesh& mesh, bool static_mesh);
  void write(Mesh& mesh, chkp_outstream& out);
  void write(std::vector<Function *> func, chkp_outstream& out);
  void write(std::vector<Vector *> vec, chkp_outstream& out);

  enum CheckpointState
  {
    CHECKPOINT, RESTART
  };
  enum RestartState
  {
    OPEN, MESH, FUNC, VEC
  };

  CheckpointState state_;
  RestartState restart_state_;

#ifdef ENABLE_MPIIO
  MPI_Offset byte_offset_;
  MPI_File in_;
#else
  std::ifstream in_;
#endif

  typedef struct
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
    uint num_shared;
#ifdef ENABLE_MPIIO
    uint offsets[5];
    uint disp[5];
#endif
  } chkp_mesh_hdr;

  chkp_mesh_hdr hdr_;

  uint n_;
  uint id_;
  real t_;
  bool hdr_initialized_;
  bool disp_initialized_;
};

}
#endif
