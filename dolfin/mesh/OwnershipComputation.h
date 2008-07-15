// Copyright (C) 2008 Niclas Jansson. 
// Licensed under the GNU LGPL Version 2.1. 

#ifndef __OWNERSHIP_COMPUTATION_H
#define __OWNERSHIP_COMPUTATION_H

#include <set>

namespace dolfin
{
  class Mesh;
  class BoundaryMesh;
  class OwnershipComputation
  {
  public:
    static void generate_ownership(Mesh& mesh);

  private:    
    static void init_edge_ownership(Mesh& mesh, BoundaryMesh& local_boundary);

    static void init_face_ownership(Mesh& mesh, BoundaryMesh& local_boundary);

    // An edge contains a pair of vertices
    typedef std::pair<uint, uint> EdgeKey;

    // A face contains a set edges
    typedef std::set<EdgeKey> FaceKey;

    // Construct a key from edge vertices
    static EdgeKey edge_key(uint id1, uint id2);
  };
}

#endif
