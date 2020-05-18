// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
//

#ifndef __DOLFIN_DOF_MAP_SET_H
#define __DOLFIN_DOF_MAP_SET_H

#include <dolfin/fem/DofMap.h>

#include <dolfin/common/types.h>

#include <ufc.h>

#include <string>
#include <vector>

namespace dolfin
{

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
  void Check(ufc::form const& form, Mesh& mesh);

  // Release dof maps
  void ReleaseAll();

  // Array of dof maps for current form
  std::vector<DofMap*> dof_map_set;

  // Mesh
  Mesh const& mesh_;

};

}

#endif

