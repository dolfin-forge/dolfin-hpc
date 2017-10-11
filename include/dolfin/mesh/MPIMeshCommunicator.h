// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm, 2007.
// Modified by Niclas Jansson, 2008.
// Modified by Balthasar Reuter, 2013.
// Modified by Aurelien Larcher, 2015-2016.
//
// The original implementations were split to ease readability, especially since
// the distribution by vertices actually does not distribute any mesh.
// The implementation was simplified and some loops were saved.
//
// First added:  2007-05-30
// Last changed: 2013-03-22

#ifndef __DOLFIN_MPI_MESH_COMMUNICATOR_H
#define __DOLFIN_MPI_MESH_COMMUNICATOR_H

#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

/// The class facilitates the transfer of a mesh between processes using MPI

struct MPIMeshCommunicator
{

  /// Distribute mesh according to a vertex-based distribution
  static void distribute(MeshValues<uint, Vertex>& distribution);

  /// Distribute mesh according to a cell-based distribution
  static void distribute(MeshValues<uint, Cell>& distribution);

};

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_MESH_COMMUNICATOR_H */
