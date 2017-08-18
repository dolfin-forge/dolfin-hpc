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

#include <dolfin/fem/DirichletBC.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
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

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, Mesh& mesh,
                         SubDomain const& sub_domain, BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain),
  g_(g),
  method_(method),
  entities_(NULL)
{
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Coefficient& g, Mesh& mesh,
                         SubDomain const& sub_domain,
                         SubSystem const& sub_system, BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain, sub_system),
  g_(g),
  method_(method),
  entities_(NULL)
{
}
//-----------------------------------------------------------------------------
DirichletBC::~DirichletBC()
{
  delete entities_;
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b,
                        BilinearForm const& form)
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
  if(form.trial_space() != form.test_space())
  {
    error("DirichletBC is implemented only for identical test and trial space");
  }

  if(entities_ == NULL || this->invalid_mesh())
  {
    if ((method_ == topological) || (method_ == geometric))
    {
      delete entities_;
      entities_ = new Array<uint>();

      // Build set of boundary facets
      uint const tdim = mesh().topology().dim();
      for (FacetIterator f(mesh()); !f.end(); ++f)
      {
        bool const on_boundary = (f->num_entities(tdim) == 1) && !f->is_shared();
        if (!this->sub_domain().inside(*f, on_boundary))
        {
          continue;
        }
        // Get cell to which facet belongs (there may be two, but pick first)
        Cell cell(mesh(), f->entities(tdim)[0]);
        entities_->push_back(cell.index());
        entities_->push_back(cell.index(*f));
      }
    }

    this->update_mesh_dependency();
  }

  // Simple check
  form.check(A,b);

  // Check compatibility of function g and the test (sub)space
  FiniteElementSpace const& space = form.trial_space();
  ufc::finite_element * fe = space.element().create_sub_element(this->sub_system());
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
void DirichletBC::computeBCTopological(_map<uint, real>& boundary_values,
                                       FiniteElementSpace const& space,
                                       SubSystem const& sub_system)
{
  // Special case
  if (entities_->size() == 0)
  {
    if(!space.mesh().is_distributed())
    {
      warning("Found no facets matching domain for boundary condition.");
    }
    return;
  }

  // Iterate over facets
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  for (Array<uint>::const_iterator it = entities_->begin();
       it != entities_->end();)
  {
    // Get cell number and local facet number
    uint const c_index = *(it++);
    uint const f_index = *(it++);

    // Create cell
    Cell cell(mesh(), c_index);
    scratch.cell.update(cell);

    // Tabulate dofs on cell for the full space dofmap
    dof_map.tabulate_dofs(cell_dofs, scratch.cell, cell);

    // Interpolate function on cell
    g_.interpolate(scratch.coefficients, scratch.cell, *scratch.finite_element,
                   cell, f_index);

    // Tabulate which dofs of the subdofmap are on the facet
    scratch.dof_map->tabulate_facet_dofs(scratch.facet_dofs, f_index);

    // Pick values for facet
    for (uint i = 0; i < scratch.dof_map->num_facet_dofs(); ++i)
    {
      uint const dof = cell_dofs[scratch.offset + scratch.facet_dofs[i]];
      real const value = scratch.coefficients[scratch.facet_dofs[i]];
      boundary_values[dof] = value;
    }
  }
  delete [] cell_dofs;
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCGeometric(_map<uint, real>& boundary_values,
                                     FiniteElementSpace const& space,
                                     SubSystem const& sub_system)
{
  // Special case
  if (entities_->size() == 0)
  {
    if(!space.mesh().is_distributed())
    {
      warning("Found no facets matching domain for boundary condition.");
    }
    return;
  }

  // Initialize facets, needed for geometric search
  uint const gdim = mesh().geometry().dim();
  uint const facet_dim = mesh().type().facet_dim();

  // Iterate over facets
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  CellType * facet_type = mesh().type().create(mesh().type().facetType());
  Point xdof;
  for (Array<uint>::const_iterator it = entities_->begin();
       it != entities_->end();)
  {
    // Get cell number and local facet number
    uint const c_index = *(it++);
    uint const f_index = *(it++);

    // Create facet
    Cell cell(mesh(), c_index);
    Facet facet(mesh(), cell.entities(facet_dim)[f_index]);

    // Loop the vertices associated with the facet
    for (VertexIterator vertex(facet); !vertex.end(); ++vertex)
    {
      // Loop the cells associated with the vertex
      for (CellIterator c(*vertex); !c.end(); ++c)
      {
        scratch.cell.update(*c);
        bool interpolated = false;

        // Tabulate dofs on cell for the full space dofmap
        dof_map.tabulate_dofs(cell_dofs, scratch.cell, *c);

        // Tabulate coordinates of dofs on cell
        scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);

        // Loop over all dofs on cell
        for (uint i = 0; i < scratch.local_dimension; ++i)
        {
          // Copy coordinates to node Point
          xdof.set(&scratch.coordinates[i][0], gdim);
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
  DofMap const& dof_map = space.dofmap();
  uint * cell_dofs = new uint[dof_map.local_dimension()];
  ScratchSpace scratch(space, sub_system);
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    scratch.cell.update(*cell);

    // Tabulate dofs on cell
    dof_map.tabulate_dofs(cell_dofs, scratch.cell, *cell);

    // Tabulate coordinates of dofs on cell
    scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);

    // Interpolate function only once and only on cells where necessary
    bool interpolated = false;

    // Loop all dofs on cell
    for (uint i = 0; i < scratch.local_dimension; ++i)
    {
      // Check if the coordinates are part of the sub domain
      if (!this->sub_domain().inside(scratch.coordinates[i], true) )
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
      uint const dof = cell_dofs[scratch.offset + i];
      real const value = scratch.coefficients[i];
      boundary_values[dof] = value;
    }
  }
  delete [] cell_dofs;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
