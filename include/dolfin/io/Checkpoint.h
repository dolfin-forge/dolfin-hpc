// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-08
// Last changed: 2010-06-20

#ifndef __CHECKPOINT_H
#define __CHECKPOINT_H

#include <fstream>
#include <string>
#include <vector>

#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/function/Function.h>

#ifdef ENABLE_MPIIo
#include <mpi.h>
#endif

namespace dolfin
{  
  class Checkpoint 
  {
  public:
    
    Checkpoint();
    ~Checkpoint();
    

    void write(std::string fname,
	       uint id, real t, Mesh& mesh,
	       std::vector<Function *> func,
	       std::vector<Vector *> vec,
	       bool static_mesh = false);

    void restart(std::string fname);

    void load(Mesh& mesh);
    void load(std::vector<Function *> func);
    void load(std::vector<Vector *> vec);

    inline bool restart() {return state == RESTART;};

    inline dolfin::uint id()
    { if(state != RESTART)
	error("Shut her down, Scotty, she's sucking mud again!");
	return _id;
    };

    inline dolfin::real restart_time()
    { if(state != RESTART) 
	error("Shut her down, Scotty, she's sucking mud again!");
      return _t;
    };

    inline void reset() 
    { state = CHECKPOINT; restart_state = OPEN; };

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




    enum CheckpointState {CHECKPOINT, RESTART};
    enum RestartState {OPEN, MESH, FUNC, VEC};

    CheckpointState state;
    RestartState restart_state;
    
#ifdef ENABLE_MPIIO
    MPI_File in;
    MPI_Offset byte_offset;
#else
    std::ifstream in;
#endif
    

    typedef struct {
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
      /*
      uint coord_offset;
      uint cells_offset;
      uint mapping_offset;
      uint ghosts_offset;
      uint shared_offset;
      */
    } chkp_mesh_hdr;


    chkp_mesh_hdr hdr;


    uint n;
    uint _id;
    real _t;
    bool hdr_initialized;
    bool disp_initialized;

  };
}
#endif
