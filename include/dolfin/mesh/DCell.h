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

  /// Dynamic cell entity to be used with DMesh
  ///
  class DCell
  {
  public:
    DCell();

    /// Check if a certain edge (defined by endpoints) is part of the cell
    bool has_edge(DVertex* v1, DVertex* v2);

    /// Local index of cell
    int id;

    /// Index of parent cell
    int parent_id;

    /// List of vertices spaning the cell
    std::vector<DVertex *> vertices;

    /// Marker for deletion
    bool deleted;
    
    /// reference number to be used for identification in bisect
    int nref;
  };
}

#endif