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
DofMapSet::DofMapSet(Form const& form, Mesh& mesh) :
    cache_(DofMapCache::instance()),
    mesh_(mesh)
{
  update(form, mesh);
}

//-----------------------------------------------------------------------------
DofMapSet::~DofMapSet()
{
  // Release all dof maps in the cache
  for (std::vector<DofMap*>::iterator it = dof_map_set.begin();
      it != dof_map_set.end(); ++it)
  {
    cache_.release_dofmap(**it);
  }
}

//-----------------------------------------------------------------------------
void DofMapSet::update(Form const& form, Mesh& mesh)
{
  // Consistency checking
#ifdef DEBUG
  Check(form, mesh);
#endif

  // Resize array of dof maps
  const uint num_arguments = form.form().rank() + form.form().num_coefficients();
  dof_map_set.resize(num_arguments);

  // Create dof maps and reuse previously computed dof maps
  for (uint i = 0; i < num_arguments; ++i)
  {
    //
    dof_map_set[i] = &cache_.acquire_dofmap(mesh, form, i);
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
  dolfin_assert(i < dof_map_set.size());
  return *dof_map_set[i];
}

//-----------------------------------------------------------------------------
void DofMapSet::Check(const ufc::form& form, Mesh& mesh)
{
  // Check that the form matches the mesh
  if (form.rank() + form.num_coefficients() > 0)
  {
    ufc::dof_map * dofmap = form.create_dof_map(0);
    if (dofmap->geometric_dimension() != mesh.geometry().dim())
    {
      error("Geometric dimension mismatch between mesh and form.");
    }
    delete dofmap;
  }
}

//-----------------------------------------------------------------------------

}


