// Copyright (C) 2007-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/SparsityPatternBuilder.h>

#include <dolfin/common/timing.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

namespace SparsityPatternBuilder
{

//-----------------------------------------------------------------------------
void build( GenericSparsityPattern& sparsity_pattern, Mesh& mesh,
            UFC& ufc, DofMapSet const& dof_map_set )
{
  message(1, "SparsityPatternBuilder: build");
  tic();

  // Initialise sparsity pattern
  sparsity_pattern.init( ufc.form.rank(),
                         ufc.global_dimensions.data(),
                         ufc.local_sizes.data() );

  // Only build for rank >= 2 (matrices and higher order tensors)
  if (ufc.form.rank() < 2)
  {
    tocd(1);
    return;
  }

  // JANPACK doesn't need any sparsity pattern information
  if( dolfin_get<std::string>("linear algebra backend") == "JANPACK" )
  {
    tocd(1);
    return;
  }

  // Build sparsity pattern for cell integrals
  if (ufc.form.has_cell_integrals())
  {
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      // Update to current cell
      ufc.cell.update(*cell);

      // Tabulate dofs for each dimension
      for (size_t i = 0; i < ufc.form.rank(); ++i)
      {
        dof_map_set[i].tabulate_dofs(ufc.dofs[i], ufc.cell);
      }

      // Fill sparsity pattern.
      sparsity_pattern.insert( ufc.local_dimensions.data(), ufc.dofs.data() );
    }
  }

  // Note: no need to iterate over exterior facets since those dofs
  // are included when tabulating dofs on all cells

  // Build sparsity pattern for interior facet integrals
  if (ufc.form.has_interior_facet_integrals())
  {
    size_t const tdim = mesh.topology_dimension();

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
      for (size_t i = 0; i < ufc.form.rank(); ++i)
      {
        const size_t offset = dof_map_set[i]().num_element_dofs();
        dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i], ufc.cell0);
        dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i] + offset, ufc.cell1);
      }

      // Fill sparsity pattern.
      sparsity_pattern.insert( ufc.macro_local_dimensions.data(),
                               ufc.macro_dofs.data() );
    }
  }

  // Build sparsity pattern for periodic facets
  // Only for square systems
  if(mesh.has_periodic_constraint())
  {
    bool has_facet_dofs = (dof_map_set[0]().num_facet_dofs() > 0);
    bool is_square = true;
    for (size_t i = 1; i < ufc.form.rank(); ++i)
    {
      has_facet_dofs |= (dof_map_set[i]().num_facet_dofs() > 0);
      is_square &= (dof_map_set[i] == dof_map_set[i - 1]);
    }

    if(has_facet_dofs)
    {
      // FIXME is this the correct space?!
      FiniteElementSpace space( mesh, ufc.finite_elements[0], dof_map_set[0]() );
      PeriodicDofsMapping const& pdm = dof_map_set[0].periodic_mapping( space );
      size_t local_dim[2];
      local_dim[0] = 1;
      local_dim[1] = pdm.max_local_dimension();
      size_t ** dofs = new size_t*[2];
      dofs[0] = new size_t[1];
      dofs[1] = new size_t[pdm.max_local_dimension()];
      for (size_t i = 0; i < pdm.num_Gdofs(); ++i)
      {
        pdm.tabulate_dofs(i, dofs[0], dofs[1], local_dim[1]);

        // Fill sparsity pattern.
        sparsity_pattern.insert(local_dim, dofs);
      }
      delete [] dofs[1];
      delete [] dofs[0];
      delete [] dofs;
    }
  }

  // Finalize sparsity pattern
  sparsity_pattern.apply();

  tocd(1);
}
//-----------------------------------------------------------------------------

} // end namespace SparsityPatternBuilder

} /* namespace dolfin */
