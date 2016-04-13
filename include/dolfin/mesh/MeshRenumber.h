// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2016.
//
// This class was fully rewritten as topology renumbering and eliminate the mesh
// layer as much as possible due to a flawed original abstraction.
// Most of the code dealing with entity distribution was consolidated in the
// from of DistributedData class.
//
// This reimplementation fixes several issues:
// - computation of boundary meshes for the distribution of faces.
// - invalid distribution of edges and faces due to invalid logical assumptions.
// - duplication of code for edges/faces.
// - inconsistency of the adjacency data.
//
// Fixes to the original implementation for the correctness of the edges/faces
// distribution were written 2013.
//

#ifndef __DOLFIN_MESH_RENUMBER_H
#define __DOLFIN_MESH_RENUMBER_H

#include <dolfin/common/types.h>

namespace dolfin
{

class MeshTopology;

/**
 *  @class  MeshRenumber
 *
 *  @brief  Provides algorithms to renumber the mesh topology:
 *          - vertices and cells exist and are only renumbered with indices
 *            contiguous per rank i.e. within the process range:
 *              [offset, offset + range size [
 *          - edges and facets need to be assigned to a rank and then numbered
 *            contiguously per rank.
 *
 */

class MeshRenumber
{

public:

  /// Renumber all mesh entities
  static bool renumber(MeshTopology& topology);

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_RENUMBER_H */
