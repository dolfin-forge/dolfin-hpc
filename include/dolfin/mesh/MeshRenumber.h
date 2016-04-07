// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
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
