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

  class DVertex
  {
  public:
    DVertex();

    int id;
    int glb_id;

    std::list<DCell *> cells;
    Point p;

    bool on_boundary;
    bool shared;
    bool ghosted;

    int owner;
  };
}

#endif