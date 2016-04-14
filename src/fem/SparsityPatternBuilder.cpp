// Copyright (C) 2007-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2007.
// Modified by Anders Logg, 2008.
//
// First added:  2007-05-24
// Last changed: 2008-02-15

#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/fem/UFC.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void SparsityPatternBuilder::build(GenericSparsityPattern& sparsity_pattern,
                                   Mesh& mesh, UFC& ufc,
                                   DofMapSet const& dof_map_set)
{
  // Initialise sparsity pattern
  bool is_distributed = mesh.is_distributed();

  if(is_distributed)
  {
    sparsity_pattern.pinit(ufc.form.rank(), ufc.global_dimensions);
    sparsity_pattern.initRange(dof_map_set[0].local_size());
  }
  else
  {
    sparsity_pattern.init(ufc.form.rank(), ufc.global_dimensions);
  }

  // Only build for rank >= 2 (matrices and higher order tensors)
  if (ufc.form.rank() < 2)
  {
    return;
  }

  // JANPACK doesn't need any sparsity pattern information
  std::string la_backend = dolfin_get("linear algebra backend");
  if(la_backend == "JANPACK")
  {
    return;
  }

  // Build sparsity pattern for cell integrals
  if (ufc.form.num_cell_integrals() != 0)
  {
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      // Update to current cell
      ufc.cell.update(*cell);

      // Tabulate dofs for each dimension
      for (uint i = 0; i < ufc.form.rank(); ++i)
      {
        dof_map_set[i].tabulate_dofs(ufc.dofs[i], ufc.cell, *cell);
      }

      // Fill sparsity pattern.
      if(is_distributed)
      {
        sparsity_pattern.pinsert(ufc.local_dimensions, ufc.dofs);
      }
      else
      {
        sparsity_pattern.insert(ufc.local_dimensions, ufc.dofs);
      }
    }
  }

  // Note: no need to iterate over exterior facets since those dofs
  // are included when tabulating dofs on all cells

  // Build sparsity pattern for interior facet integrals
  if (ufc.form.num_interior_facet_integrals() != 0)
  {
    uint const tdim = mesh.topology().dim();

    // Compute facets and facet - cell connectivity if not already computed
    mesh.init(tdim - 1);
    mesh.init(tdim - 1, tdim);

    for (FacetIterator facet(mesh); !facet.end(); ++facet)
    {
      // Check if we have an interior facet
      if (facet->num_entities(tdim) != 2)
      {
        continue;
      }

      // Get cells incident with facet
      Cell cell0(mesh, facet->entities(tdim)[0]);
      Cell cell1(mesh, facet->entities(tdim)[1]);

      // Update to current pair of cells
      ufc.cell0.update(cell0);
      ufc.cell1.update(cell1);

      // Tabulate dofs for each dimension on macro element
      for (uint i = 0; i < ufc.form.rank(); ++i)
      {
        const uint offset = dof_map_set[i].local_dimension();
        dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i], ufc.cell0, cell0);
        dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i] + offset, ufc.cell1, cell1);
      }

      // Fill sparsity pattern.
      if(is_distributed)
      {
        sparsity_pattern.pinsert(ufc.macro_local_dimensions, ufc.macro_dofs);
      }
      else
      {
        sparsity_pattern.insert(ufc.macro_local_dimensions, ufc.macro_dofs);
      }
    }
  }

  // Build sparsity pattern for periodic facets
  // Only for square systems
  if(mesh.has_periodic_constraint())
  {
    bool has_facet_dofs = (dof_map_set[0].num_facet_dofs() > 0);
    bool is_square = true;
    for (uint i = 1; i < ufc.form.rank(); ++i)
    {
      has_facet_dofs |= (dof_map_set[i].num_facet_dofs() > 0);
      is_square &= (dof_map_set[i] == dof_map_set[i - 1]);
    }

    if(has_facet_dofs)
    {
      PeriodicDofsMapping const& pdm = dof_map_set[0].periodic_mapping();
      uint local_dim[2];
      local_dim[0] = 1;
      local_dim[1] = pdm.max_local_dimension();
      uint ** dofs = new uint*[2];
      dofs[0] = new uint[1];
      dofs[1] = new uint[pdm.max_local_dimension()];
      for (uint i = 0; i < pdm.num_Gdofs(); ++i)
      {
        pdm.tabulate_dofs(i, dofs[0], dofs[1], local_dim[1]);

        // Fill sparsity pattern.
        if(is_distributed)
        {
          sparsity_pattern.pinsert(local_dim, dofs);
        }
        else
        {
          sparsity_pattern.insert(local_dim, dofs);
        }
      }
      delete [] dofs[1];
      delete [] dofs[0];
      delete [] dofs;
    }
  }

  // Finalize sparsity pattern
  sparsity_pattern.apply();
}
//-----------------------------------------------------------------------------

}

