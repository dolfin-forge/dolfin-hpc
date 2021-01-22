// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/timing.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/PeriodicDofsMapping.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/UFCHalo.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Scalar.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/log/log.h>
#include <dolfin/main/OpenMP.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>

#include <memory>
#include <vector>

namespace dolfin
{

namespace Assembler
{

/// Assemble tensor from given (UFC) form, coefficients and sub domains.
/// This is the main assembly function in DOLFIN. All other assembly functions
/// end up calling this function.
///
/// The MeshFunction arguments can be used to specify assembly over subdomains
/// of the mesh cells, exterior facets and interior facets.
/// Either a null pointer or an empty MeshFunction may be used to specify that
/// the tensor should be assembled over the entire set of cells or facets.
void assemble(GenericTensor& A, const Form& form,
              std::vector<Coefficient*> const& coefficients,
              DofMapSet const& dofmaps,
              MeshValues<size_t, Cell> const* cell_domains,
              MeshValues<size_t, Facet> const* exterior_facet_domains,
              MeshValues<size_t, Facet> const* interior_facet_domains,
              bool reset_tensor = true);

// Assemble over cells
void assembleCells(GenericTensor& A, std::vector<Coefficient*> const& coefficients,
                   DofMapSet const& dofmaps, UFC& data,
                   MeshValues<size_t, Cell> const* domains);

// Assemble over exterior facets
void assembleExteriorFacets(GenericTensor& A,
                            std::vector<Coefficient*> const& coefficients,
                            DofMapSet const& dofmaps, UFC& data,
                            MeshValues<size_t, Facet> const* domains);

// Assemble over interior facets
void assembleInteriorFacets(GenericTensor& A,
                            std::vector<Coefficient*> const& coefficients,
                            DofMapSet const& dofmaps, UFC& data,
                            MeshValues<size_t, Facet> const* domains);

// Bogus-assemble periodic contributions
void initializePeriodicDofs(GenericTensor& A,
                            std::vector<Coefficient*> const& coefficients,
                            DofMapSet const& dofmaps, UFC& data,
                            MeshValues<size_t, Facet> const* domains);

// Initialize global tensor
void initGlobalTensor(GenericTensor& A, DofMapSet const& dofmaps, UFC& ufc,
                      bool reset_tensor);

//-----------------------------------------------------------------------------
void assemble( GenericTensor & A, Form & form, bool reset_tensor )
{
  OPENMP_PRAGMA( parallel )
  assemble( A, form, form.coefficients(), form.dofmaps(),
            nullptr, nullptr, nullptr, reset_tensor );
}
//-----------------------------------------------------------------------------
void assemble(GenericTensor& A, Form& form,
                         SubDomain const& sub_domain, bool reset_tensor)
{
  Mesh& mesh = form.mesh();

  // Extract cell domains
  MeshValues<size_t, Cell>* cell_domains = nullptr;

  // Extract facet domains
  MeshValues<size_t, Facet>* facet_domains = nullptr;

OPENMP_PRAGMA( master )
  {
    if ( form.has_cell_integrals() )
    {
      cell_domains = new MeshValues<size_t, Cell>(mesh);
      (*cell_domains) = 1;
      sub_domain.mark(*cell_domains, 0);
    }

    if ( form.has_exterior_facet_integrals() or
         form.has_interior_facet_integrals() )
    {
      facet_domains = new MeshValues<size_t, Facet>(mesh);
      (*facet_domains) = 1;
      sub_domain.mark(*facet_domains, 0);
    }
  }

  // Assemble
  assemble(A, form, form.coefficients(), form.dofmaps(),
           cell_domains, facet_domains, facet_domains, reset_tensor);

  // Delete domains
OPENMP_PRAGMA( master )
  {
    delete cell_domains;
    delete facet_domains;
  }
}
//-----------------------------------------------------------------------------
void assemble(GenericTensor& A, Form& form,
                         MeshValues<size_t, Cell> const& cell_domains,
                         MeshValues<size_t, Facet> const& exterior_facet_domains,
                         MeshValues<size_t, Facet> const& interior_facet_domains,
                         bool reset_tensor)
{
  assemble(A, form, form.coefficients(), form.dofmaps(), &cell_domains,
           &exterior_facet_domains, &interior_facet_domains, reset_tensor);
}
//-----------------------------------------------------------------------------
void assemble(GenericTensor& A, const Form& form,
                         std::vector<Coefficient*> const& coefficients,
                         DofMapSet const& dofmaps,
                         MeshValues<size_t, Cell> const* cell_domains,
                         MeshValues<size_t, Facet> const* exterior_facet_domains,
                         MeshValues<size_t, Facet> const* interior_facet_domains,
                         bool reset_tensor)
{
  // Check arguments
OPENMP_PRAGMA( master )
  {
    if(reset_tensor)
    {
      form.check(coefficients);
    }
  }

  // Create data structure for local assembly data
  UFC ufc(form);

  // Initialize global tensor
OPENMP_PRAGMA( master )
  {
    initGlobalTensor(A, dofmaps, ufc, reset_tensor);


    // Update all ghost degrees of freedom
    for ( Coefficient * coeff : coefficients )
    {
      coeff->sync();
    }
  }
OPENMP_PRAGMA( flush )
OPENMP_PRAGMA( barrier )

  // Assemble over cells
  assembleCells(A, coefficients, dofmaps, ufc, cell_domains);

  // Assemble over exterior facets
  assembleExteriorFacets(A, coefficients, dofmaps, ufc, exterior_facet_domains);

  // Assemble over interior facets
  assembleInteriorFacets(A, coefficients, dofmaps, ufc, interior_facet_domains);

  // Bogus-assemble periodic dofs
  initializePeriodicDofs(A, coefficients, dofmaps, ufc, exterior_facet_domains);

  // Finalise assembly of global tensor
OPENMP_PRAGMA( master )
  A.apply();
OPENMP_PRAGMA( barrier )
}
//-----------------------------------------------------------------------------
void assembleCells(GenericTensor& A,
                              std::vector<Coefficient*> const& coefficients,
                              DofMapSet const& dofmaps,
                              UFC& ufc,
                              MeshValues<size_t, Cell> const* domains)
{
  if ( not ufc.form.has_cell_integrals() )
  {
    return;
  }

  message(1,"Assembler: cells");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  size_t const N = mesh.num_cells();
  size_t const form_rank = ufc.form.rank();
  size_t const coef_size = coefficients.size();
  ufc::cell_integral * integral = ufc.cell_integrals[0];

  CellIterator it(mesh);
OPENMP_PRAGMA( for )
  for (size_t i = 0; i < N; ++i)
  {
    Cell& cell = it[i];

    // Get integral for sub domain (if any)
    if ((domains != nullptr) && domains->size() > 0)
    {
      size_t const domain = (*domains)(cell);
      // FIXME is this correct?!
      if (domain <= ufc.form.max_cell_subdomain_id())
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
    for (size_t c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], cell);
    }

    // Tabulate dofs for each dimension
    for (size_t d = 0; d < form_rank; ++d)
    {
      dofmaps[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate cell tensor
    // FIXME last argument (cell_orientation) needs an actual value
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell.coordinates.data(), 0);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void assembleExteriorFacets(GenericTensor& A,
                                       std::vector<Coefficient*> const& coefficients,
                                       DofMapSet const& dofmaps,
                                       UFC& ufc,
                                       MeshValues<size_t, Facet> const* domains)
{
  if ( not ufc.form.has_exterior_facet_integrals() )
  {
    return;
  }

  message(1,"Assembler: exterior facets");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  size_t const tdim = mesh.topology_dimension();
  BoundaryMesh& exterior_boundary = mesh.exterior_boundary();
  size_t const N = exterior_boundary.num_cells();
  if (N == 0)
  {
    return;
  }
  size_t const form_rank = ufc.form.rank();
  size_t const coef_size = coefficients.size();
  ufc::exterior_facet_integral * integral = ufc.exterior_facet_integrals[0];

  FacetIterator it(mesh);
  CellIterator  c0(mesh);
OPENMP_PRAGMA( for )
  for (size_t i = 0; i < N; ++i)
  {
    // Get mesh facet corresponding to boundary cell
    Facet& facet = it[exterior_boundary.facet_index(i)];

    // Get integral for sub domain (if any)
    if ((domains != nullptr) && domains->size() > 0)
    {
      size_t const domain = (*domains)(facet);
      // FIXME is this correct?!
      if (domain <= ufc.form.max_exterior_facet_subdomain_id())
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
    size_t const local_facet = cell.index(facet);

    // Update to current cell
    ufc.cell.update(cell);

    // Interpolate coefficients on cell
    for (size_t c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], cell, local_facet);
    }

    // Tabulate dofs for each dimension
    for (size_t d = 0; d < form_rank; ++d)
    {
      dofmaps[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate exterior facet tensor
    // FIXME last argument (cell_orientation) needs an actual value
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell.coordinates.data(),
                              local_facet, 0);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void assembleInteriorFacets(GenericTensor& A,
                                       std::vector<Coefficient*> const& coefficients,
                                       DofMapSet const& dofmaps,
                                       UFC& ufc,
                                       MeshValues<size_t, Facet> const* domains)
{
  if ( not ufc.form.has_interior_facet_integrals() )
  {
    return;
  }

  message(1,"Assembler: interior facets");
  tic();

  Mesh& mesh = dofmaps[0].mesh();
  size_t const tdim = mesh.topology_dimension();
  size_t const N = mesh.size(mesh.type().facet_dim());
  size_t const form_rank = ufc.form.rank();
  size_t const coef_size = coefficients.size();
  ufc::interior_facet_integral * integral = ufc.interior_facet_integrals[0];

  // Halo data structure caching macro element coefficients and dofs
  UFCHalo halo(ufc, coefficients, dofmaps);

  FacetIterator it(mesh);
  CellIterator  c0(mesh);
  CellIterator  c1(mesh);
OPENMP_PRAGMA( for )
  for (size_t i = 0; i < N; ++i)
  {
    Facet& facet = it[i];

    // Get integral for sub domain (if any)
    if ((domains != nullptr) && domains->size() > 0)
    {
      size_t const domain = (*domains)(facet);
      // FIXME is this correct?!
      if (domain <= ufc.form.max_interior_facet_subdomain_id())
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
      for (size_t c = 0; c < coef_size; ++c)
      {
        coefficients[c]->interpolate(ufc.macro_w[c], ufc.cell0, *ufc.coefficient_elements[c], cell0, ufc.facet0);
        size_t const offset = ufc.coefficient_elements[c]->space_dimension();
        coefficients[c]->interpolate(ufc.macro_w[c] + offset, ufc.cell1, *ufc.coefficient_elements[c], cell1, ufc.facet1);
      }

      // Tabulate dofs for each dimension on cell1
      for (size_t d = 0; d < form_rank; ++d)
      {
        dofmaps[d].tabulate_dofs(ufc.macro_dofs[d], ufc.cell0);
        size_t const offset = ufc.local_dimensions[d];
        dofmaps[d].tabulate_dofs(ufc.macro_dofs[d] + offset, ufc.cell1);
      }

      // FIXME last two arguments (cell_orientation) need an actual value
      integral->tabulate_tensor(ufc.macro_A, ufc.macro_w,
                                ufc.cell0.coordinates.data(),
                                ufc.cell1.coordinates.data(),
                                ufc.facet0, ufc.facet1, 0, 0);

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

      // FIXME last two arguments (cell_orientation) need an actual value
      integral->tabulate_tensor(ufc.macro_A, halo.macro_w,
                                halo.cell0.coordinates.data(),
                                halo.cell1.coordinates.data(),
                                halo.facet0, halo.facet1, 0, 0);

      // Add entries to global tensor
      A.add(ufc.macro_A, ufc.macro_local_dimensions, ufc.macro_dofs);
    }
  }

  tocd(1);
}
//-----------------------------------------------------------------------------
void initializePeriodicDofs(GenericTensor& A,
                            std::vector<Coefficient*> const&,
                            DofMapSet const& dofmaps,
                            UFC& ufc,
                            MeshValues<size_t, Facet> const*)
{
  if(!dofmaps[0].mesh().has_periodic_constraint())
  {
    return;
  }

  // Add zero at periodic dofs to allocate entries
  /// @todo This could be fixed by a modification of the assembler's behaviour
  if(A.rank() == 2)
  {
    //
    Matrix& matA = static_cast<Matrix&>(A);
    // FIXME is this the correct space?!
    Mesh & mesh = const_cast<Mesh&>( dofmaps.mesh() );
    FiniteElementSpace space( mesh, ufc.finite_elements[0], dofmaps[0] );
    PeriodicDofsMapping const& pdm = dofmaps[0].periodic_mapping( space );
    real * block = new real[pdm.max_local_dimension() + 1];
    std::fill_n(block, pdm.max_local_dimension() + 1, 0.0);
    size_t irow = 0;
    size_t * jcols = new size_t[pdm.max_local_dimension() + 1];
    std::fill_n(jcols, pdm.max_local_dimension() + 1, 0.0);
    size_t ncols = 0;
    for (size_t i = 0; i < pdm.num_Gdofs(); ++i)
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
void initGlobalTensor(GenericTensor& A, DofMapSet const& dofmaps,
                      UFC& ufc, bool reset_tensor)
{
  if (A.rank() == 0)
  {
    if ( reset_tensor )
    {
      A.zero();
    }

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

} // namespace Assembler

} // namespace dolfin
