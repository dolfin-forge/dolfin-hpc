// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008
// Modified by Ola Skavhaug, 2007
// Modified by Niclas Jansson, 2008-2010.
//
// First added:  2007-01-17
// Last changed: 2010-03-18

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Scalar.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/UFCHalo.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/common/timing.h>

#include <map>
#include <memory>

namespace dolfin
{

//-----------------------------------------------------------------------------
Assembler::Assembler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Assembler::~Assembler()
{
}
//-----------------------------------------------------------------------------
void Assembler::assemble(GenericTensor& A, Form& form, bool reset_tensor)
{
#pragma omp parallel
  assemble(A, form, form.coefficients(), form.dofmaps(), 0, 0, 0, reset_tensor);
}
//-----------------------------------------------------------------------------
void Assembler::assemble(GenericTensor& A, Form& form,
                         SubDomain const& sub_domain, bool reset_tensor)
{
  Mesh& mesh = form.mesh();
  uint const tdim = mesh.topology().dim();

  // Extract cell domains
  MeshFunction<uint>* cell_domains = NULL;

  // Extract facet domains
  MeshFunction<uint>* facet_domains = NULL;

#pragma omp master
  {
    if (form.num_cell_integrals() > 0)
    {
      cell_domains = new MeshFunction<uint>(mesh, tdim);
      (*cell_domains) = 1;
      sub_domain.mark(*cell_domains, 0);
    }

    if (form.num_exterior_facet_integrals() > 0 ||
        form.num_interior_facet_integrals() > 0)
    {
      facet_domains = new MeshFunction<uint>(mesh, mesh.type().facet_dim());
      (*facet_domains) = 1;
      sub_domain.mark(*facet_domains, 0);
    }
  }

  // Assemble
  assemble(A, form, form.coefficients(), form.dofmaps(),
           cell_domains, facet_domains, facet_domains, reset_tensor);

  // Delete domains
#pragma omp master
  {
    delete cell_domains;
    delete facet_domains;
  }
}
//-----------------------------------------------------------------------------
void Assembler::assemble(GenericTensor& A, Form& form,
                         MeshFunction<uint> const& cell_domains,
                         MeshFunction<uint> const& exterior_facet_domains,
                         MeshFunction<uint> const& interior_facet_domains,
                         bool reset_tensor)
{
  assemble(A, form, form.coefficients(), form.dofmaps(), &cell_domains,
           &exterior_facet_domains, &interior_facet_domains, reset_tensor);
}
//-----------------------------------------------------------------------------
void Assembler::assemble(GenericTensor& A, const Form& form,
                         Array<Coefficient*> const& coefficients,
                         DofMapSet const& dofmaps,
                         MeshFunction<uint> const* cell_domains,
                         MeshFunction<uint> const* exterior_facet_domains,
                         MeshFunction<uint> const* interior_facet_domains,
                         bool reset_tensor)
{
  // Check arguments
#pragma omp master
  {
    if(reset_tensor)
    {
      form.check(coefficients);
    }
  }

  // Create data structure for local assembly data
  UFC ufc(form);

  // Initialize global tensor
#pragma omp master
  {
    initGlobalTensor(A, dofmaps, ufc, reset_tensor);


    // Update all ghost degrees of freedom
    for (uint i = 0; i < coefficients.size(); ++i)
    {
      coefficients[i]->sync();
    }
  }
#pragma omp flush
#pragma omp barrier

  // Assemble over cells
  assembleCells(A, coefficients, dofmaps, ufc, cell_domains);

  // Assemble over exterior facets
  assembleExteriorFacets(A, coefficients, dofmaps, ufc, exterior_facet_domains);

  // Assemble over interior facets
  assembleInteriorFacets(A, coefficients, dofmaps, ufc, interior_facet_domains);

  // Bogus-assemble periodic dofs
  initializePeriodicDofs(A, coefficients, dofmaps, ufc, exterior_facet_domains);

  // Finalise assembly of global tensor
#pragma omp master
  A.apply();
#pragma omp barrier
}
//-----------------------------------------------------------------------------
void Assembler::assembleCells(GenericTensor& A,
                              Array<Coefficient*> const& coefficients,
                              DofMapSet const& dofmaps,
                              UFC& ufc,
                              MeshFunction<uint> const* domains) const
{
  if (ufc.form.num_cell_integrals() == 0)
  {
    return;
  }

  message(1,"Assembler : cells");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  uint const N = mesh.num_cells();
  uint const form_rank = ufc.form.rank();
  uint const coef_size = coefficients.size();
  ufc::cell_integral * integral = ufc.cell_integrals[0];

  CellIterator it(mesh);
#pragma omp for
  for (uint i = 0; i < N; ++i)
  {
    Cell& cell = it[i];

    // Get integral for sub domain (if any)
    if ((domains != NULL) && domains->size() > 0)
    {
      uint const domain = (*domains)(cell);
      if (domain < ufc.form.num_cell_integrals())
      {
          integral = ufc.cell_integrals[domain];
      }
      else
      {
        continue;
      }
    }

    // Update to current cell
    ufc.cell.update(cell);

    // Interpolate coefficients on cell
    for (uint c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], cell);
    }

    // Tabulate dofs for each dimension
    for (uint d = 0; d < form_rank; ++d)
    {
      dofmaps[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate cell tensor
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void Assembler::assembleExteriorFacets(GenericTensor& A,
                                       Array<Coefficient*> const& coefficients,
                                       DofMapSet const& dofmaps,
                                       UFC& ufc,
                                       MeshFunction<uint> const* domains) const
{
  if (ufc.form.num_exterior_facet_integrals() == 0)
  {
    return;
  }

  message(1,"Assembler : exterior facets");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  uint const tdim = mesh.topology().dim();
  BoundaryMesh& exterior_boundary = mesh.exterior_boundary();
  uint const N = exterior_boundary.num_cells();
  if (N == 0)
  {
    return;
  }
  uint const form_rank = ufc.form.rank();
  uint const coef_size = coefficients.size();
  ufc::exterior_facet_integral * integral = ufc.exterior_facet_integrals[0];

  FacetIterator it(mesh);
  CellIterator  c0(mesh);
#pragma omp for
  for (uint i = 0; i < N; ++i)
  {
    // Get mesh facet corresponding to boundary cell
    Facet& facet = it[exterior_boundary.facet_index(i)];

    // Get integral for sub domain (if any)
    if ((domains != NULL) && domains->size() > 0)
    {
      uint const domain = (*domains)(facet);
      if (domain < ufc.form.num_exterior_facet_integrals())
      {
        integral = ufc.exterior_facet_integrals[domain];
      }
      else
      {
        continue;
      }
    }

    // Get mesh cell to which mesh facet belongs
    Cell& cell = c0[facet.entities(tdim)[0]];

    // Get local index of facet with respect to the cell
    uint const local_facet = cell.index(facet);

    // Update to current cell
    ufc.cell.update(cell);

    // Interpolate coefficients on cell
    for (uint c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], cell, local_facet);
    }

    // Tabulate dofs for each dimension
    for (uint d = 0; d < form_rank; ++d)
    {
      dofmaps[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate exterior facet tensor
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell, local_facet);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void Assembler::assembleInteriorFacets(GenericTensor& A,
                                       Array<Coefficient*> const& coefficients,
                                       DofMapSet const& dofmaps,
                                       UFC& ufc,
                                       MeshFunction<uint> const* domains) const
{
  if (ufc.form.num_interior_facet_integrals() == 0)
  {
    return;
  }

  message(1,"Assembler : interior facets");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  uint const tdim = mesh.topology().dim();
  uint const N = mesh.size(mesh.type().facet_dim());
  uint const form_rank = ufc.form.rank();
  uint const coef_size = coefficients.size();
  ufc::interior_facet_integral * integral = ufc.interior_facet_integrals[0];

  // Halo data structure caching macro element coefficients and dofs
  UFCHalo halo(ufc, coefficients, dofmaps);

  FacetIterator it(mesh);
  CellIterator  c0(mesh);
  CellIterator  c1(mesh);
#pragma omp for
  for (uint i = 0; i < N; ++i)
  {
    Facet& facet = it[i];

    // Get integral for sub domain (if any)
    if ((domains != NULL) && domains->size() > 0)
    {
      uint const domain = (*domains)(facet);
      if (domain < ufc.form.num_interior_facet_integrals())
      {
        integral = ufc.interior_facet_integrals[domain];
      }
      else
      {
        continue;
      }
    }

    // Check if we have a local interior facet
    if (facet.num_entities(tdim) == 2)
    {
      // Get cells incident with facet, local index of facet
      Cell& cell0 = c0[facet.entities(tdim)[0]];
      ufc.facet0 = cell0.index(facet);
      ufc.cell0.update(cell0);

      Cell& cell1 = c1[facet.entities(tdim)[1]];
      ufc.facet1 = cell1.index(facet);
      ufc.cell1.update(cell1);

      // Interpolate coefficients on cell1
      for (uint c = 0; c < coef_size; ++c)
      {
        coefficients[c]->interpolate(ufc.macro_w[c], ufc.cell0, *ufc.coefficient_elements[c], cell0, ufc.facet0);
        uint const offset = ufc.coefficient_elements[c]->space_dimension();
        coefficients[c]->interpolate(ufc.macro_w[c] + offset, ufc.cell1, *ufc.coefficient_elements[c], cell1, ufc.facet1);
      }

      // Tabulate dofs for each dimension on cell1
      for (uint d = 0; d < form_rank; ++d)
      {
        dofmaps[d].tabulate_dofs(ufc.macro_dofs[d], ufc.cell0);
        uint const offset = ufc.local_dimensions[d];
        dofmaps[d].tabulate_dofs(ufc.macro_dofs[d] + offset, ufc.cell1);
      }

      integral->tabulate_tensor(ufc.macro_A, ufc.macro_w, ufc.cell0, ufc.cell1,
                                ufc.facet0, ufc.facet1);

      // Add entries to global tensor
      A.add(ufc.macro_A, ufc.macro_local_dimensions, ufc.macro_dofs);
    }
    // Interprocess facet
    else if (facet.is_shared())
    {
      // Contributions from cell0 are restored from the halo data while
      // contributions from cell1 are fetched from adjacent ranks.
      // Implementation updates pointers to coordinates, but has to copy dofs
      // and coefficients until data structures are reworked.
      halo.update(facet);

      integral->tabulate_tensor(ufc.macro_A, halo.macro_w, halo.cell0,
                                halo.cell1, halo.facet0, halo.facet1);

      // Add entries to global tensor
      A.add(ufc.macro_A, ufc.macro_local_dimensions, ufc.macro_dofs);
    }
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void Assembler::initializePeriodicDofs(GenericTensor& A,
                                       Array<Coefficient*> const& coefficients,
                                       DofMapSet const& dofmaps,
                                       UFC& data,
                                       MeshFunction<uint> const* domains) const
{
  if(!dofmaps[0].mesh().has_periodic_constraint())
  {
    return;
  }

  // Add zero at periodic dofs to allocate entries
  //FIXME: This could be fixed by a modification of the assembler's behaviour
  if(A.rank() == 2)
  {
    //
    Matrix& matA = static_cast<Matrix&>(A);
    PeriodicDofsMapping const& pdm = dofmaps[0].periodic_mapping();
    real * block = new real[pdm.max_local_dimension() + 1];
    std::fill_n(block, pdm.max_local_dimension() + 1, 0.0);
    uint irow = 0;
    uint * jcols = new uint[pdm.max_local_dimension() + 1];
    std::fill_n(jcols, pdm.max_local_dimension() + 1, 0.0);
    uint ncols = 0;
    for (uint i = 0; i < pdm.num_Gdofs(); ++i)
    {
      pdm.tabulate_dofs(i, &irow, jcols, ncols);
      jcols[ncols] = irow;
      matA.add(block, 1, &irow, ncols, jcols);
    }
    delete[] jcols;
    delete[] block;
  }
}
//-----------------------------------------------------------------------------
void Assembler::initGlobalTensor(GenericTensor& A, DofMapSet const& dofmaps,
                                 UFC& ufc, bool reset_tensor) const
{
  if (A.rank() == 0)
  {
    return;
  }

  if (reset_tensor || A.size(0) == 0)
  {
    GenericSparsityPattern * sparsity_pattern = A.factory().createPattern();
    SparsityPatternBuilder::build(*sparsity_pattern, dofmaps[0].mesh(), ufc,
                                  dofmaps);
    A.init(*sparsity_pattern);
    delete sparsity_pattern;
  }
  else
  {
    if((A.rank() > 0 && A.size(0) != ufc.global_dimensions[0]) ||
       (A.rank() > 1 && A.size(1) != ufc.global_dimensions[1]))
    {
      error("Assembler : dimensions of linear system do not match spaces.");
    }
    A.zero();
  }

}
//-----------------------------------------------------------------------------

}
