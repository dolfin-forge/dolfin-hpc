// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_RENUMBER_H
#define __DOLFIN_MESH_RENUMBER_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;

/**
 *  @class  MeshRenumber
 *
 *  @brief  Provides algorithms to renumber mesh entities
 *  .
 */

class MeshRenumber
{

public:

  /// Renumber all mesh entities
  static bool renumber(Mesh& mesh);

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_RENUMBER_H */
