// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Aurélien Larcher, 2013. (implement dofmap cache)
//
// First added:  2007-01-17
// Last changed: 2007-05-24

#ifndef __DOF_MAP_SET_H
#define __DOF_MAP_SET_H

#include <map>
#include <vector>
#include <string>
#include <ufc.h>

#include <dolfin/common/types.h>
#include "DofMap.h"

namespace dolfin
{

class DofMapCache;
class Form;
class Mesh;
class UFC;

/// This class provides storage and caching of (precomputed) dof
/// maps and enables reuse of already computed dof maps with equal
/// signatures.

class DofMapSet
{
public:

  /// Create set of dof maps
  DofMapSet(Form const& form, Mesh& mesh);

  /// Destructor
  ~DofMapSet();

  /// Update set of dof maps for given form
  void update(Form const& form, Mesh& mesh);

  /// Return number of dof maps
  uint size() const;

  /// Return the mesh on which dof maps have been initialized
  Mesh const& mesh() const;

  /// Return dof map for argument function i
  DofMap& operator[](uint i) const;

private:

  // Consistency checking
  void Check(const ufc::form& form, Mesh& mesh);

  // Array of dof maps for current form
  std::vector<DofMap*> dof_map_set;

  // Global cache of precomputed dof maps
  DofMapCache& cache_;

  // Mesh
  Mesh const& mesh_;

};

}

#endif

