// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

// The original implementations were split to ease readability, especially since
// the distribution by vertices actually does not distribute any mesh.
// The implementation was simplified and some loops were saved.

#ifndef __DOLFIN_MPI_MESH_COMMUNICATOR_H
#define __DOLFIN_MPI_MESH_COMMUNICATOR_H

#include <dolfin/mesh/MeshData.h>

namespace dolfin
{

/// The namesapce facilitates the transfer of a mesh between processes using MPI

namespace MPIMeshCommunicator
{

/// Distribute mesh according to a vertex-based distribution
void distribute( MeshValues< uint, Vertex > & dist );

/// Distribute mesh according to a cell-based distribution
void distribute( MeshValues< uint, Cell > & dist, MeshData * D = nullptr );

/// Check mesh entity distribution
template < class E >
void check( Mesh & Mesh );

/// Check mesh distribution
void check( Mesh & Mesh );

} /* namespace MPIMeshCommunicator */

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_MESH_COMMUNICATOR_H */
