// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#ifndef __D_CELL_H
#define __D_CELL_H

#include <vector>

namespace dolfin
{
  class DVertex;

  /// Representation of a cell within a dynamic mesh DMesh
  
  class DCell
  {
  public:
    DCell();
    bool has_edge(DVertex* v1, DVertex* v2);

    int id;
    int parent_id;

    std::vector<DVertex *> vertices;

    bool deleted;
    
    int nref;
  };
}

#endif