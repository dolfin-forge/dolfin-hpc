// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2016.
//
// This class was first fixed in 2014 due to an invalid assumption causing
// incorrect computation of the global facet map then modified to support
// hypercube meshes.
// In 2015 it was cleaned up and rewritten to reduce communication.
//
// First added:  2008-02-25
// Last changed: 2008-02-25

#ifndef __DOLFIN_GLOBAL_FACET_MAP_H
#define __DOLFIN_GLOBAL_FACET_MAP_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Facet;
class Mesh;

class GlobalFacetMap
{

public:

  ///
  GlobalFacetMap(Mesh& mesh);

  ///
  ~GlobalFacetMap();

  ///
  bool is_global(Facet& facet);

  ///
  bool is_shared(Facet& facet);

  ///
  void disp() const;

private:

  ///
  void init();

  //
  Mesh& mesh_;

  //
  uint const tdim_;

  // Shared facets on local boundary
  _set<uint> shared_facets_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_GLOBAL_FACET_MAP_H */
