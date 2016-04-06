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
 *  @brief  Provides algorithms to renumber the mesh topology.
 *  .
 */

class MeshRenumber
{

public:

  /// Renumber all mesh entities
  static bool renumber(MeshTopology& topology);

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_RENUMBER_H */
