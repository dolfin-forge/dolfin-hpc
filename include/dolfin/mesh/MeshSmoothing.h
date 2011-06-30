// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Jeannette Spuhler, Rodrigo Vilela De Abreu and Kaspar Muller 2011.
// First added:  2008-07-16
// Last changed: 2011-06-30

#ifndef __MESH_SMOOTHING_H
#define __MESH_SMOOTHING_H

#include <dolfin/mesh/MeshSmoothData.h>
#include <map>

namespace dolfin
{
  
  class Mesh;

  /// This class implements mesh smoothing. The coordinates of
  /// internal vertices are updated by local averaging.

  class MeshSmoothing
  {

  public:

    static void smooth(Mesh& mesh);

    static void smooth(Mesh& mesh, MeshSmoothData& smooth_data);
    
  private:

    static void smooth_common(Mesh& mesh, MeshSmoothData& smooth_data);

  };
}

#endif
