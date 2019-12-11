// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __MESH_SMOOTHING_H
#define __MESH_SMOOTHING_H

namespace dolfin
{

  class Mesh;

  /// This class implements mesh smoothing. The coordinates of
  /// internal vertices are updated by local averaging.

namespace MeshSmoothing
{

void smooth( Mesh & mesh );

} // namespace MeshSmoothing

} // namespace dolfin

#endif
