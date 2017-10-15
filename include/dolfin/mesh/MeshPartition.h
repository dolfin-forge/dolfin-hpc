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
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Mesh;

/// This class provides a set of functions to partition a Mesh

class MeshPartition
{

public:

  /// Partition a mesh into pe_size partitions in parallel
  static void partition(MeshValues<uint, Cell>& partitions);

  /// Partition a mesh into pe_size partitions in parallel with
  /// weights on vertices ^H^H^H on the *fucking* *cells*
  static void partition(MeshValues<uint, Cell>& partitions, MeshValues<uint, Cell>& weight);

  /// Partition a mesh based on coordinates
  static void partition_geom(MeshValues<uint, Vertex>& partitions);
  
};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_PARTITIONING_H */
