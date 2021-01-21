// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
//

#ifndef __DOLFIN_DOF_MAP_SET_H
#define __DOLFIN_DOF_MAP_SET_H

#include <dolfin/common/types.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/ufc/ufc.h>

#include <string>

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
  auto size() const -> uint;

  /// Return the mesh on which dof maps have been initialized
  auto mesh() const -> Mesh const&;

  /// Return dof map for argument function i
  auto operator[](uint i) const -> DofMap&;

private:

  // Consistency checking
  void Check(ufc::form const& form, Mesh& mesh);

  // Release dof maps
  void ReleaseAll();

  // Array of dof maps for current form
  Array<DofMap*> dof_map_set;

  // Mesh
  Mesh const& mesh_;

};

//-----------------------------------------------------------------------------
inline auto DofMapSet::size() const -> uint
{
  return dof_map_set.size();
}

//-----------------------------------------------------------------------------
inline auto DofMapSet::mesh() const -> Mesh const &
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline auto DofMapSet::operator[]( uint i ) const -> DofMap &
{
  dolfin_assert( dof_map_set.size() > 0 );
  dolfin_assert( i < dof_map_set.size() );
  return *dof_map_set[i];
}
}

#endif

