// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by: Magnus Vikstrøm, 2007.
// Modified by: Niclas Jansson, 2008.
// Modified by: Balthasar Reuter, 2013.
//
// First added:  2007-05-30
// Last changed: 2008-01-14

#ifndef __MPI_MESH_COMMUNICATOR_H
#define __MPI_MESH_COMMUNICATOR_H

#include "MeshFunction.h"

namespace dolfin
{
  class Mesh;
  class MPI;

  /// The class facilitates the transfer of a mesh between processes using MPI
  
  class MPIMeshCommunicator
  {
  public:
    
    /// Constructor
    MPIMeshCommunicator();

    /// Destructor
    ~MPIMeshCommunicator();

    /// Distribute mesh according to a mesh function
    static void distribute(Mesh& mesh, MeshFunction<uint>& distribution);
    
    /// Distribute mesh according to mesh function and preserve cell markers
    static void distribute(Mesh& mesh, MeshFunction<uint>& distribution, 
			   MeshFunction<bool>& old_cell_marker,
			   MeshFunction<bool>& cell_marker);

    /// Distribute mesh according to mesh function and preserve mesh functions
    static void distribute(Mesh& mesh, MeshFunction<uint>& distribution,
         Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
          >& cell_functions);
    
  private:
    static void distributeCommon(Mesh& mesh, MeshFunction<uint>& distribution, 
				 MeshFunction<bool>* old_cell_marker,
				 MeshFunction<bool>* cell_marker);

    static void distributeCommon(Mesh& mesh, MeshFunction<uint>& distribution,
         Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
          >& cell_functions);
  };
}

#endif
