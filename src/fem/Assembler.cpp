// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008
// Modified by Ola Skavhaug, 2007
// Modified by Niclas Jansson, 2008-2010.
//
// First added:  2007-01-17
// Last changed: 2010-03-18

#include <memory>
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
#include <dolfin/mesh/MeshData.h>
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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
Assembler::Assembler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Assembler::Assembler(Mesh& mesh)
{
  // Do nothing
  warning("Assembler(Mesh& mesh) is deprecated, use Assembler().");
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
                         const SubDomain& sub_domain, bool reset_tensor)
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
      facet_domains = new MeshFunction<uint>(mesh, tdim - 1);
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
dolfin::real Assembler::assemble(Form& form)
{
  Scalar value;
  assemble(value, form, true);
  return value;
}
//-----------------------------------------------------------------------------
dolfin::real Assembler::assemble(Form& form,
                                 const SubDomain& sub_domain)
{
  Scalar value;
  assemble(value, form, sub_domain, true);
  return value;
}
//-----------------------------------------------------------------------------
dolfin::real Assembler::assemble(Form& form,
                                 MeshFunction<uint> const& cell_domains,
                                 MeshFunction<uint> const& exterior_facet_domains,
                                 MeshFunction<uint> const& interior_facet_domains)
{
  Scalar value;
  assemble(value, form,
           cell_domains, exterior_facet_domains, interior_facet_domains, true);
  return value;
}
//-----------------------------------------------------------------------------
void Assembler::assemble(GenericTensor& A, const Form& form,
                         Array<Coefficient*> const& coefficients,
                         const DofMapSet& dof_map_set,
                         const MeshFunction<uint>* cell_domains,
                         const MeshFunction<uint>* exterior_facet_domains,
                         const MeshFunction<uint>* interior_facet_domains,
                         bool reset_tensor)
{
  // Check arguments
#pragma omp master
  {
    if(reset_tensor)
    {
      form.check_coefficients(coefficients);
    }
  }

  // Create data structure for local assembly data
  UFC ufc(form, form.mesh(), dof_map_set);

  // Initialize global tensor
#pragma omp master
  {
    initGlobalTensor(A, dof_map_set, ufc, reset_tensor);


    // Update all ghost points
    for (uint i = 0; i < coefficients.size(); ++i)
    {
      coefficients[i]->sync();
    }
  }
#pragma omp flush
#pragma omp barrier

  // Assemble over cells
  assembleCells(A, coefficients, dof_map_set, ufc, cell_domains);

  // Assemble over exterior facets
  assembleExteriorFacets(A, coefficients, dof_map_set, ufc, exterior_facet_domains);

  // Assemble over interior facets
  assembleInteriorFacets(A, coefficients, dof_map_set, ufc, interior_facet_domains);

  // Bogus-assemble periodic dofs
  initializePeriodicDofs(A, coefficients, dof_map_set, ufc, exterior_facet_domains);

  // Finalise assembly of global tensor
#pragma omp master
  A.apply();
#pragma omp barrier
}
//-----------------------------------------------------------------------------
void Assembler::assembleCells(GenericTensor& A,
                              Array<Coefficient*> const& coefficients,
                              const DofMapSet& dof_map_set,
                              UFC& ufc,
                              const MeshFunction<uint>* domains) const
{
  // Skip assembly if there are no cell integrals
  if (ufc.form.num_cell_integrals() == 0)
    return;

  Mesh& mesh = dof_map_set[0].mesh();

  // Cell integral
  ufc::cell_integral* integral = ufc.cell_integrals[0];

  // Assemble over cells
#ifndef NO_PROGRESS_BAR
  Progress p(progressMessage(A.rank(), "cells"), mesh.numCells());
#endif
  //  for (CellIterator cell(mesh); !cell.end(); ++cell)
#pragma omp for
  uint const num_cells = mesh.numCells();
  uint const coef_size = coefficients.size();
  uint const form_rank = ufc.form.rank();
  MeshDistributedData& distdata = mesh.distdata();
  for (uint i = 0; i < num_cells; ++i)
  {
    Cell cell(mesh, i);

    // Get integral for sub domain (if any)
    if (domains && domains->size() > 0)
    {
      const uint domain = (*domains)(cell);
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
    ufc.update(cell, distdata);

    // Interpolate coefficients on cell
    for (uint c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], cell);
    }

    // Tabulate dofs for each dimension
    for (uint d = 0; d < form_rank; ++d)
    {
      dof_map_set[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate cell tensor
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);

#ifndef NO_PROGRESS_BAR
    p++;
#endif
  }

}
//-----------------------------------------------------------------------------
void Assembler::assembleExteriorFacets(GenericTensor& A,
                                       Array<Coefficient*> const& coefficients,
                                       const DofMapSet& dof_map_set,
                                       UFC& ufc,
                                       const MeshFunction<uint>* domains) const
{
  // Skip assembly if there are no exterior facet integrals
  if (ufc.form.num_exterior_facet_integrals() == 0)
    return;

  Mesh& mesh = dof_map_set[0].mesh();
  uint const tdim = mesh.topology().dim();

  // Exterior facet integral
  ufc::exterior_facet_integral* integral = ufc.exterior_facet_integrals[0];

  BoundaryMesh& exterior_boundary = mesh.exterior_boundary();
  if(exterior_boundary.numCells()  == 0) return;
  MeshFunction<uint>* cell_map = exterior_boundary.data().meshFunction("cell map");

  dolfin_assert(cell_map);

  // Assemble over exterior facets (the cells of the boundary)
#ifndef NO_PROGRESS_BAR
  dolfin_assert(exterior_boundary.numCells());
  Progress p(progressMessage(A.rank(), "exterior facets"), exterior_boundary.numCells());
#endif
  //  for (CellIterator boundary_cell(*boundary); !boundary_cell.end(); ++boundary_cell)
#pragma omp for
  uint const ext_num_cells = exterior_boundary.numCells();
  uint const coef_size = coefficients.size();
  uint const form_rank = ufc.form.rank();
  MeshDistributedData& distdata = mesh.distdata();
  for (uint i = 0; i < ext_num_cells; i++)
  {
    // Get mesh facet corresponding to boundary cell
    Facet mesh_facet(mesh, (*cell_map).get(i));

    // Get integral for sub domain (if any)
    if (domains && domains->size() > 0)
    {
      const uint domain = (*domains)(mesh_facet);
      if (domain < ufc.form.num_exterior_facet_integrals())
      {
        integral = ufc.exterior_facet_integrals[domain];
      }
      else
      {
        continue;
      }
    }

    // Get mesh cell to which mesh facet belongs (pick first, there is only one)
    dolfin_assert(mesh_facet.num_entities(tdim) == 1);
    Cell mesh_cell(mesh, mesh_facet.entities(tdim)[0]);

    // Get local index of facet with respect to the cell
    const uint local_facet = mesh_cell.index(mesh_facet);

    // Update to current cell
    ufc.update(mesh_cell, distdata);

    // Interpolate coefficients on cell
    for (uint c = 0; c < coef_size; ++c)
    {
      coefficients[c]->interpolate(ufc.w[c], ufc.cell, *ufc.coefficient_elements[c], mesh_cell, local_facet);
    }

    // Tabulate dofs for each dimension
    for (uint d = 0; d < form_rank; ++d)
    {
      dof_map_set[d].tabulate_dofs(ufc.dofs[d], ufc.cell);
    }

    // Tabulate exterior facet tensor
    integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell, local_facet);

    // Add entries to global tensor
    A.add(ufc.A, ufc.local_dimensions, ufc.dofs);

#ifndef NO_PROGRESS_BAR
    p++;
#endif

  }
}
//-----------------------------------------------------------------------------
void Assembler::assembleInteriorFacets(GenericTensor& A,
                                       Array<Coefficient*> const& coefficients,
                                       const DofMapSet& dof_map_set,
                                       UFC& ufc,
                                       const MeshFunction<uint>* domains) const
{
  // Skip assembly if there are no interior facet integrals
  if (ufc.form.num_interior_facet_integrals() == 0)
  {
    return;
  }

  Mesh& mesh = dof_map_set[0].mesh();
  uint const tdim = mesh.topology().dim();

  // Interior facet integral
  ufc::interior_facet_integral* integral = ufc.interior_facet_integrals[0];

  // Compute facets and facet - cell connectivity if not already computed
  uint const facet_dim = tdim - 1;
  mesh.init(facet_dim);
  mesh.init(facet_dim, tdim);
  mesh.order();

  // Assemble over interior facets (the facets of the mesh)
#ifndef NO_PROGRESS_BAR
  dolfin_assert(mesh.numFacets());
  Progress p(progressMessage(A.rank(), "interior facets"), mesh.numFacets());
#endif

  // The halo data structure caches macro element coefficients and dofs for each
  // shared facet
  UFCHalo halo(ufc, coefficients, dof_map_set);

  //
  for (FacetIterator facet(mesh); !facet.end(); ++facet)
  {
    // Check if we have an interior facet
    if (facet->num_entities(tdim) != 2)
    {
#ifndef NO_PROGRESS_BAR
      p++;
#endif
      continue;
    }

    // Get integral for sub domain (if any)
    if (domains && domains->size() > 0)
    {
      const uint domain = (*domains)(*facet);
      if (domain < ufc.form.num_interior_facet_integrals())
      {
        integral = ufc.interior_facet_integrals[domain];
      }
      else
      {
        continue;
      }
    }

    if(!facet->is_shared())
    {
      // Get cells incident with facet, local index of facet
      Cell cell0(mesh, facet->entities(tdim)[0]);
      ufc.facet0 = cell0.index(*facet);
      ufc.cell0.update(cell0);

      // Contributions from cell1 are local
      Cell cell1(mesh, facet->entities(tdim)[1]);
      ufc.facet1 = cell1.index(*facet);
      ufc.cell1.update(cell1);

      // Interpolate coefficients on cell1
      for (uint c = 0; c < coefficients.size(); ++c)
      {
        coefficients[c]->interpolate(ufc.macro_w[c], ufc.cell0, *ufc.coefficient_elements[c], cell0, ufc.facet0);
        uint const offset = ufc.coefficient_elements[c]->space_dimension();
        coefficients[c]->interpolate(ufc.macro_w[c] + offset, ufc.cell1, *ufc.coefficient_elements[c], cell1, ufc.facet1);
      }

      // Tabulate dofs for each dimension on cell1
      for (uint d = 0; d < ufc.form.rank(); ++d)
      {
        dof_map_set[d].tabulate_dofs(ufc.macro_dofs[d], ufc.cell0, cell0.index());
        uint const offset = ufc.local_dimensions[d];
        dof_map_set[d].tabulate_dofs(ufc.macro_dofs[d] + offset, ufc.cell1, cell1.index());
      }

      integral->tabulate_tensor(ufc.macro_A, ufc.macro_w, ufc.cell0, ufc.cell1,
                                ufc.facet0, ufc.facet1);
    }
    else
    {
      // Contributions from cell0 is restored from data stored in the halo data
      // while contribution from cell1 are fetched from adjacent ranks
      // Implementation updates pointers to coordinates, but has to copy dofs
      // and coefficients until data structures are reworked.
      halo.update(*facet);

      integral->tabulate_tensor(ufc.macro_A, halo.macro_w, halo.cell0,
                                halo.cell1, halo.facet0, halo.facet1);

    }

    // Add entries to global tensor
    A.add(ufc.macro_A, ufc.macro_local_dimensions, ufc.macro_dofs);

#ifndef NO_PROGRESS_BAR
    p++;
#endif
  }

}
//-----------------------------------------------------------------------------
void Assembler::initializePeriodicDofs(GenericTensor& A,
                                       const Array<Coefficient*>& coefficients,
                                       const DofMapSet& dof_map_set,
                                       UFC& data,
                                       const MeshFunction<uint>* domains) const
{
  if(!dof_map_set[0].mesh().has_periodic_constraint())
  {
    return;
  }

  // Add zero at periodic dofs to allocate entries
  //FIXME: This could be fixed by a modification of the assembler's behaviour
  if(A.rank() == 2)
  {
    //
    Matrix& matA = static_cast<Matrix&>(A);
    PeriodicDofsMapping const& pdm = dof_map_set[0].periodic_mapping();
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
void Assembler::initGlobalTensor(GenericTensor& A, const DofMapSet& dof_map_set,
                                 UFC& ufc, bool reset_tensor) const
{

  if (reset_tensor || A.size(0) == 0)
  {
    GenericSparsityPattern * sparsity_pattern = A.factory().createPattern();
    SparsityPatternBuilder::build(*sparsity_pattern, dof_map_set[0].mesh(), ufc,
                                  dof_map_set);
    A.init(*sparsity_pattern);
    delete sparsity_pattern;
  }
  else
  {
    if((A.rank() > 0 && A.size(0) != ufc.global_dimensions[0])
        || (A.rank() > 1 && A.size(1) != ufc.global_dimensions[1]))
    {
      error("Dimensions of linear system do not match discrete spaces.");
    }
    A.zero();
  }

}
//-----------------------------------------------------------------------------
std::string Assembler::progressMessage(uint rank, std::string integral_type) const
{
  std::stringstream s;
  s << "Assembling (rank " << MPI::processNumber()<< " ) ";

  switch (rank)
  {
  case 0:
    s << "scalar value over ";
    break;
  case 1:
    s << "vector over ";
    break;
  case 2:
    s << "matrix over ";
    break;
  default:
    s << "rank " << rank << " tensor over ";
    break;
  }

  s << integral_type;

  return s.str();
}

//-----------------------------------------------------------------------------

}
