// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __MESH_RENUMBER_H
#define __MESH_RENUMBER_H

#include "Face.h"
#include <dolfin/common/Array.h>
#include <set>

namespace dolfin
{

class Mesh;

/**
 *  @class  MeshRenumber
 *
 *  @brief  Provides algorithms to renumbers mesh entities and store the indices
 *          into one array of local map and one array of global maps.
 *          The entities numbering are stored at index given by 3D definition:
 *            - 0: vertex, filled by renumber_vertices(Mesh& mesh),
 *            - 1: edge, filled by renumber_edges(Mesh& mesh),
 *            - 2: face, filled by renumber_faces(Mesh& mesh),
 *            - 3: cell, filled by renumber_cells(Mesh& mesh).
 *          In the 2D, the same indices apply but faces are ignored.
 *          In the 1D, the same indices apply but faces are ignored.
 *          If P1 optimizations are enabled only vertices and cells are
 *          renumbered.
 */

class MeshRenumber
{
public:

  /// Renumber all mesh entities
  static void renumber(Mesh& mesh);

  /// Renumber mesh vertices
  static void renumber_vertices(Mesh& mesh);

  /// Renumber mesh edges
  static void renumber_edges(Mesh& mesh);

  /// Renumber mesh faces
  static void renumber_faces(Mesh& mesh);

  /// Renumber mesh cells
  static void renumber_cells(Mesh& mesh);

private:

  /// An edge contains a pair of vertices
  typedef std::pair<uint, uint> EdgeKey;

  /// A face contains a set edges
  typedef std::set<EdgeKey> FaceKey;

  /// Construct a key from edge vertices
  static EdgeKey edge_key(uint id1, uint id2);

  /// Construct a key from face vertices
  static FaceKey face_key(Face& f);

  static void send_buffer_face(Array<uint>& send_buff, Mesh& mesh, Face& f);

};

}

#endif
