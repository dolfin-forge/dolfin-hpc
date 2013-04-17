// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#ifndef __D_VERTEX_H
#define __D_VERTEX_H

#include <list>

#include <dolfin/mesh/Point.h>

namespace dolfin
{
  class DCell;

  /// Dynamic vertex entity to be used with DMesh
  ///
  class DVertex
  {
  public:
    DVertex();

    /// Local index of vertex
    int id;

    /// Global index of vertex
    int glb_id;

    /// List of cells containing the vertex
    std::list<DCell *> cells;

    /// Vertex coordinates as Point object
    Point p;

    /// Indicator if vertex is on process boundary
    bool on_boundary;

    /// Indicator if vertex is shared
    bool shared;

    /// Indicator if vertex is ghosted
    bool ghosted;

    /// Rank of owning process
    int owner;
  };
}

#endif