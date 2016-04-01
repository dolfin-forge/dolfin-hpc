// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Kristian Oelgaard, 2007
// Modified by Martin Sandve Alnes, 2008
// Modified by Niclas Jansson, 2008-2015
// Modified by Aurélien Larcher, 2014
//
// First added:  2007-04-10
// Last changed: 2014-04-15

#include <dolfin/common/constants.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/fem/DirichletBC.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, Mesh& mesh,
                         const SubDomain& sub_domain, BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain),
  g_(g),
  method_(method)
{
  initFromSubDomain(sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, MeshFunction<uint>& sub_domains,
                         uint sub_domain, BCMethod method) :
  BoundaryCondition("Dirichlet", sub_domains, sub_domain),
  g_(g),
  method_(method)
{
  initFromMeshFunction(sub_domains, sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, Mesh& mesh,
                         const SubDomain& sub_domain,
                         const SubSystem& sub_system, BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain, sub_system),
  g_(g),
  method_(method)
{
  initFromSubDomain(sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, MeshFunction<uint>& sub_domains,
                         uint sub_domain, const SubSystem& sub_system,
                         BCMethod method) :
  BoundaryCondition("Dirichlet", sub_domains, sub_domain, sub_system),
  g_(g),
  method_(method)
{
  initFromMeshFunction(sub_domains, sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::~DirichletBC()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form)
{
  apply(A, b, 0, form);
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b,
                        GenericVector const& x, BilinearForm const& form)
{
  apply(A, b, &x, form);
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b,
                        GenericVector const* x, BilinearForm const& form)
{
  // Get local data for application of boundary conditions
  FiniteElementSpace const& space = form.trial_space();

  if(form.trial_space() != form.test_space())
  {
    error("DirichletBC is implemented only for identical test and trial space");
  }

  if(this->invalid_mesh())
  {
    warning("Mesh topology and geometry have changed");
    facets_.clear();
    if(this->has_geometrical_sub_domain())
    {
      initFromSubDomain(this->sub_domain());
    }
    else
    {
      initFromMeshFunction(this->sub_domain_markers(),
                           this->sub_domain_index());
    }
    this->update_mesh_dependency();
    message("Facet map for Dirichlet boundary condition recomputed");
  }

  // Simple check
  form.check(A,b);

  // Check compatibility of function g and the test (sub)space
  ufc::finite_element * fe =
      space.element().create_sub_element(this->sub_system());
  if((fe->value_rank() != g_.rank())||(fe->value_dimension(0)!=g_.dim(0)))
  {
    error("Rank and/or value dimension mismatch between function and space.\n"
          "Function : rank = %d, dim = %d; Space : rank = %d, dim = %d.",
          g_.rank(), g_.dim(0), fe->value_rank(), fe->value_dimension(0));
  }
  delete fe;

  // A map to hold the mapping from boundary dofs to boundary values
  _map<uint, real> boundary_values;

  // Compute dofs and values
  switch (method_)
    {
    case topological:
      computeBCTopological(boundary_values, space, this->sub_system());
      break;
    case geometric:
      computeBCGeometric(boundary_values, space, this->sub_system());
      break;
    case pointwise:
      computeBCPointwise(boundary_values, space, this->sub_system());
      break;
    default:
      error("Unknown method for application of boundary conditions.");
      break;
    }

  // Copy boundary value data to arrays
  uint * dofs = new uint[boundary_values.size()];
  real * values = new real[boundary_values.size()];
  _map<uint, real>::const_iterator boundary_value;
  uint i = 0;
  for (boundary_value = boundary_values.begin();
       boundary_value != boundary_values.end(); ++boundary_value)
  {
    dofs[i] = boundary_value->first;
    values[i++] = boundary_value->second;
  }

  // Modify boundary values for nonlinear problems
  if (x)
  {
    real * x_values = new real[boundary_values.size()];
    x->get(x_values, boundary_values.size(), dofs);
    for (uint i = 0; i < boundary_values.size(); ++i)
    {
      values[i] -= x_values[i];
    }
    delete[] x_values;
  }

  message("Applying boundary conditions to linear system");

  // Modify RHS vector (b[i] = value)
  b.set(values, boundary_values.size(), dofs);

  // Modify linear system (A_ii = 1)
  A.ident(boundary_values.size(), dofs);

  // Clear temporary arrays
  delete[] dofs;
  delete[] values;

  // Finalise changes to A
  A.apply();

  // Finalise changes to b
  b.apply();
}
//-----------------------------------------------------------------------------
void DirichletBC::initFromSubDomain(SubDomain const& sub_domain)
{
  dolfin_assert(facets_.size() == 0);

  if(mesh().size(0) == 0)
  {
    error("Provided mesh is empty");
  }

  // FIXME: This can be made more efficient, we should be able to
  // FIXME: extract the facets without first creating a MeshFunction on
  // FIXME: the entire mesh and then extracting the subset. This is done
  // FIXME: mainly for convenience (we may reuse mark() in SubDomain).

  // Create mesh function for sub domain markers on facets
  const uint dim = mesh().topology().dim();
  mesh().init(dim - 1);
  MeshFunction<uint> sub_domains(mesh(), dim - 1);

  // Mark everything as sub domain 1
  sub_domains = 1;

  // Mark the sub domain as sub domain 0
  sub_domain.mark(sub_domains, 0);

  // Initialize from mesh function
  initFromMeshFunction(sub_domains, 0);
}
//-----------------------------------------------------------------------------
void DirichletBC::initFromMeshFunction(MeshFunction<uint> const& sub_domains,
                                       uint sub_domain)
{
  dolfin_assert(facets_.size() == 0);

  // Make sure we have the facet - cell connectivity
  uint const dim = mesh().topology().dim();
  mesh().init(dim - 1, dim);

  // Build set of boundary facets
  for (FacetIterator facet(mesh()); !facet.end(); ++facet)
  {
    // Skip facets not on this boundary
    if (sub_domains(*facet) != sub_domain)
    {
      continue;
    }

    // Get cell to which facet belongs (there may be two, but pick first)
    Cell cell(mesh(), facet->entities(dim)[0]);

    // Get local index of facet with respect to the cell
    const uint facet_number = cell.index(*facet);

    // Copy data
    facets_.push_back(std::pair<uint, uint>(cell.index(), facet_number));
  }
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCTopological(_map<uint, real>& boundary_values,
                                       FiniteElementSpace const& space,
                                       SubSystem const& sub_system)
{
  // Special case
  if (facets_.size() == 0)
  {
    if(!space.mesh().is_distributed())
    {
      warning("Found no facets matching domain for boundary condition.");
    }
    return;
  }

  // Iterate over facets
#ifndef NO_PROGRESS_BAR
  Progress p("Computing Dirichlet boundary values, topological search", facets_.size());
#endif
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  for (uint f = 0; f < facets_.size(); ++f)
  {
    // Get cell number and local facet number
    uint cell_number = facets_[f].first;
    uint facet_number = facets_[f].second;

    // Create cell
    Cell cell(mesh(), cell_number);
    scratch.cell.update(cell);

    // Tabulate dofs on cell for the full space dofmap
    dof_map.tabulate_dofs(cell_dofs, scratch.cell, cell.index());

    // Interpolate function on cell
    //FIXME: DISABLED stupid facet thing, breaks FacetNormal !
    g_.interpolate(scratch.coefficients, scratch.cell, *scratch.finite_element, cell);

    // Tabulate which dofs of the subdofmap are on the facet
    scratch.dof_map->tabulate_facet_dofs(scratch.facet_dofs, facet_number);

    // Pick values for facet
    for (uint i = 0; i < scratch.dof_map->num_facet_dofs(); ++i)
    {
      uint const dof = cell_dofs[scratch.offset + scratch.facet_dofs[i]];
      real const value = scratch.coefficients[scratch.facet_dofs[i]];
      boundary_values[dof] = value;
    }

#ifndef NO_PROGRESS_BAR
    p++;
#endif
  }
  delete [] cell_dofs;
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCGeometric(_map<uint, real>& boundary_values,
                                     FiniteElementSpace const& space,
                                     SubSystem const& sub_system)
{
  // Special case
  if (facets_.size() == 0)
  {
    if(!space.mesh().is_distributed())
    {
      warning("Found no facets matching domain for boundary condition.");
    }
    return;
  }

  // Initialize facets, needed for geometric search
  message("Computing facets, needed for geometric application of boundary conditions.");
  uint const tdim = mesh().topology().dim();
  uint const facet_dim = tdim - 1;
  mesh().init(facet_dim);

  // Iterate over facets
#ifndef NO_PROGRESS_BAR
  Progress p("Computing Dirichlet boundary values, geometric search", facets_.size());
#endif
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  CellType * facet_type = mesh().type().create(mesh().type().facetType());
  Point xdof;
  for (uint f = 0; f < facets_.size(); ++f)
  {
    // Get cell number and local facet number
    uint cell_number = facets_[f].first;
    uint facet_number = facets_[f].second;

    // Create facet
    Cell cell(mesh(), cell_number);
    Facet facet(mesh(), cell.entities(facet_dim)[facet_number]);

    // Loop the vertices associated with the facet
    for (VertexIterator vertex(facet); !vertex.end(); ++vertex)
    {
      // Loop the cells associated with the vertex
      for (CellIterator c(*vertex); !c.end(); ++c)
      {
        scratch.cell.update(*c);
        bool interpolated = false;

        // Tabulate dofs on cell for the full space dofmap
        dof_map.tabulate_dofs(cell_dofs, scratch.cell, c->index());

        // Tabulate coordinates of dofs on cell
        scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);

        // Loop over all dofs on cell
        for (uint i = 0; i < scratch.local_dimension; ++i)
        {
          // Copy coordinates to node Point
          std::memcpy(&xdof[0], &scratch.coordinates[i][0], tdim*sizeof(real));
          // Check if the coordinates are on current facet and thus on boundary
          if (!facet_type->intersects(facet, xdof))
          {
            continue;
          }

          if(!interpolated)
          {
            interpolated = true;
            // Interpolate function on cell for the given (sub)element
            g_.interpolate(scratch.coefficients, scratch.cell, *scratch.finite_element, *c);
          }

          // Set boundary value
          uint const dof = cell_dofs[scratch.offset + i];
          real const value = scratch.coefficients[i];
          boundary_values[dof] = value;
        }

      }
    }
  }
  delete facet_type;
  delete [] cell_dofs;
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCPointwise(_map<uint, real>& boundary_values,
                                     FiniteElementSpace const& space,
                                     SubSystem const& sub_system)
{
  // Iterate over cells
#ifndef NO_PROGRESS_BAR
  Progress p("Computing Dirichlet boundary values, pointwise search", mesh().numCells());
#endif
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    scratch.cell.update(*cell);

    // Tabulate dofs on cell
    dof_map.tabulate_dofs(cell_dofs, scratch.cell, cell->index());

    // Tabulate coordinates of dofs on cell
    scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);

    // Interpolate function only once and only on cells where necessary
    bool interpolated = false;

    // Loop all dofs on cell
    for (uint i = 0; i < scratch.local_dimension; ++i)
    {
      // Check if the coordinates are part of the sub domain
      if ( !this->sub_domain().inside(scratch.coordinates[i], true) )
      {
        continue;
      }

      if(!interpolated)
      {
        interpolated = true;
        // Interpolate function on cell
        g_.interpolate(scratch.coefficients, scratch.cell, *scratch.finite_element, *cell);
      }

      // Set boundary value
      const uint dof = cell_dofs[scratch.offset + i];
      const real value = scratch.coefficients[i];
      boundary_values[dof] = value;
    }

#ifndef NO_PROGRESS_BAR
    p++;
#endif
  }
  delete [] cell_dofs;
}

//-----------------------------------------------------------------------------

}
