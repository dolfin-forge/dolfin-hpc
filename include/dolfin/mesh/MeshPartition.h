// Copyright (C) 2007 Magnus Vikstrom
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_PARTITIONING_H
#define __DOLFIN_MESH_PARTITIONING_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

/// This class provides a set of functions to partition a Mesh

namespace MeshPartition
{

/// Partition a mesh into pe_size partitions in parallel
void partition( MeshValues< uint, Cell > & partitions );

/// Partition a mesh into pe_size partitions in parallel with
/// weights on vertices on the cells
void partition( MeshValues< uint, Cell > & partitions,
                MeshValues< uint, Cell > & weight );

/// Partition a mesh based on coordinates
void partition_geom( MeshValues< uint, Vertex > & partitions );

} /* namespace MeshPartition */

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_PARTITIONING_H */
