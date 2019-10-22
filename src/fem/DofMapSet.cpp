// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Martin Sandve Alnes, 2008.
// Modified by Aurélien Larcher, 2013. (implement dofmap cache)
//
// First added:  2007-01-17
// Last changed: 2008-05-08

#include <dolfin/log/dolfin_log.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/fem/Form.h>

#include <ufc.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
DofMapSet::DofMapSet(Form const&, Mesh& mesh) :
    mesh_(mesh)
{
}

//-----------------------------------------------------------------------------
DofMapSet::~DofMapSet()
{
  // Release all dof maps in the cache
  ReleaseAll();
}

//-----------------------------------------------------------------------------
void DofMapSet::update(Form const& form, Mesh& mesh)
{
  // Consistency checking
#if DEBUG
  Check(form, mesh);
#endif

  // Release previously acquired dof maps if any
  ReleaseAll();

  // Resize array of dof maps
  uint const num_arguments = form.rank() + form.num_coefficients();
  dof_map_set.resize(num_arguments);

  // Create dof maps and reuse previously computed dof maps
  for (uint i = 0; i < num_arguments; ++i)
  {
    //
    dof_map_set[i] = &(DofMap::acquire(mesh, form, i));
  }
}

//-----------------------------------------------------------------------------
uint DofMapSet::size() const
{
  return dof_map_set.size();
}

//-----------------------------------------------------------------------------
Mesh const& DofMapSet::mesh() const
{
  return mesh_;
}

//-----------------------------------------------------------------------------
DofMap& DofMapSet::operator[](uint i) const
{
  dolfin_assert(dof_map_set.size() > 0);dolfin_assert(i < dof_map_set.size());
  return *dof_map_set[i];
}

//-----------------------------------------------------------------------------
void DofMapSet::Check(ufc::form const& form, Mesh& mesh)
{
  // Check that the form matches the mesh
  if (form.rank() + form.num_coefficients() > 0)
  {
    ufc::dofmap * dofmap = form.create_dofmap(0);
    if (dofmap->geometric_dimension() != mesh.geometry_dimension())
    {
      error("Geometric dimension mismatch between mesh and form.");
    }
    delete dofmap;
  }
}

//-----------------------------------------------------------------------------
void DofMapSet::ReleaseAll()
{
  // Release all dof maps in the cache
  for (std::vector<DofMap*>::iterator it = dof_map_set.begin();
      it != dof_map_set.end(); ++it)
  {
    DofMap::release(**it);
  }
}

//-----------------------------------------------------------------------------

}

