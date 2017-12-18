// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#ifndef __DOLFIN_D_VERTEX_H
#define __DOLFIN_D_VERTEX_H

#include <dolfin/mesh/Vertex.h>

#include <list>

namespace dolfin
{

class DCell;

/// Dynamic vertex entity to be used with DMesh

struct DVertex
{
  static uint const UNDEF = DOLFIN_UINT_MAX;

  DVertex();

  DVertex(Vertex const& v);

  /// Local index of vertex
  uint id;

  /// Global index of vertex
  uint glb_id;

  /// List of cells containing the vertex
  std::list<DCell *> cells;

  /// Vertex coordinates as Point object
  Point p;

  /// Marker for deletion
  bool deleted;

  /// Indicator if vertex is shared
  bool shared;

  /// Indicator if vertex is ghosted
  bool ghosted;

  /// Rank of owning process
  uint owner;

  /// Adjacent processes for boundary vertices
  _set<uint> shared_adj;
};

} /* namespace dolfin */

#endif /* __DOLFIN_D_VERTEX_H */
