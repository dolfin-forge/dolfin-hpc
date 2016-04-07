// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by: Magnus Vikstrøm, 2007.
// Modified by: Niclas Jansson, 2008.
// Modified by: Balthasar Reuter, 2013.
//
// First added:  2007-05-30
// Last changed: 2013-03-22

#ifndef __DOLFIN_MPI_MESH_COMMUNICATOR_H
#define __DOLFIN_MPI_MESH_COMMUNICATOR_H

#include <dolfin/mesh/MeshFunction.h>

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

};

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_MESH_COMMUNICATOR_H */
