// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-02-25
// Last changed: 2008-02-25

#ifndef _GLOBAL_FACET_MAP_H
#define _GLOBAL_FACET_MAP_H

#include "Mesh.h"
#include "Facet.h"
#include <map>

namespace dolfin
{
  class GlobalFacetMap
  {
  public:
    GlobalFacetMap(Mesh& mesh);
    ~GlobalFacetMap();

    void init();
    bool globalFacet(Facet& facet);
    
  private:

    void findGlobal2D();
    void findGlobal3D();

    // Map of processor global facets, 
    // stores only the facets with shared vertices
    _map<uint, bool> global_facet;
    Mesh& _mesh;
  };
}
#endif
