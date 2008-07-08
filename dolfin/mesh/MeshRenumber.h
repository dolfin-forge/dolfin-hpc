// Copyright (C) 2008 Niclas Jansson. 
// Licensed under the GNU LGPL Version 2.1. 



#ifndef __MESH_RENUMBER_H
#define __MESH_RENUMBER_H

#include "Mesh.h"


namespace dolfin
{
  class MeshRenumber
  {
  public:

    static void renumber(Mesh& mesh);

    static void renumber_vertices(Mesh& mesh);
    
    static void renumber_edges(Mesh& mesh);
    
    static void renumber_faces(Mesh& mesh);
    
    static void renumber_cells(Mesh& mesh);
    
  private:

    // An edge contains a pair of vertices
    typedef std::pair<uint, uint> EdgeKey;

    // A face contains a set edges
    typedef std::set<EdgeKey> FaceKey;

    // Construct a key from edge vertices
    static EdgeKey edge_key(uint id1, uint id2);
    

  };

}

#endif
