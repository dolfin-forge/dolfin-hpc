// Copyright (C) 2007 Magnus Vikstrom
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2008-2015.
//
// First added:  2007-04-24
// Last changed: 2015-01-05

#ifndef __DOLFIN_MESH_PARTITIONING_H
#define __DOLFIN_MESH_PARTITIONING_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;
template<class T> class MeshFunction;

/// This class provides a set of functions to partition a Mesh

class MeshPartition
{

public:

  /// Partition a mesh into pe_size partitions in parallel
  static void partition(Mesh& mesh, MeshFunction<uint>& partitions);

  /// Partition a mesh into pe_size partitions in parallel with
  /// weights on vertices
  static void partition(Mesh& mesh, MeshFunction<uint>& partitions,
                        MeshFunction<uint>& weight);

  /// Partition a mesh based on coordinates
  static void partition_geom(Mesh& mesh, MeshFunction<uint>& partitions);
  
};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_PARTITIONING_H */
