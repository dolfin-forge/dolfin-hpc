// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Kristian Oelgaard, 2007
// Modified by Martin Sandve Alnes, 2008
// Modified by Niclas Jansson, 2008
//
// First added:  2007-04-10
// Last changed: 2008-07-07

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
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/fem/DirichletBC.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Function& g, Mesh& mesh, const SubDomain& sub_domain,
                         BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain),
  g_(g),
  method_(method)
{
  initFromSubDomain(sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Function& g, MeshFunction<uint>& sub_domains,
                         uint sub_domain, BCMethod method) :
  BoundaryCondition("Dirichlet", sub_domains, sub_domain),
  g_(g),
  method_(method)
{
  initFromMeshFunction(sub_domains, sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Function& g, Mesh& mesh, const SubDomain& sub_domain,
                         const SubSystem& sub_system, BCMethod method) :
  BoundaryCondition("Dirichlet", mesh, sub_domain, sub_system),
  g_(g),
  method_(method)
{
  initFromSubDomain(sub_domain);
}
//-----------------------------------------------------------------------------
DirichletBC::DirichletBC(Function& g, MeshFunction<uint>& sub_domains,
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
void DirichletBC::apply(GenericMatrix& A, GenericVector& b, Form const& form)
{
  apply(A, b, 0, form);
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b,
                        GenericVector const& x, Form const& form)
{
  apply(A, b, &x, form);
}
//-----------------------------------------------------------------------------
void DirichletBC::apply(GenericMatrix& A, GenericVector& b,
                        GenericVector const* x, Form const& form)
{
  // Get local data for application of boundary conditions
  BoundaryCondition::LocalData& data = this->updateLocalData(form);

  // Simple check
  const uint N = data.dof_map->global_dimension();
  if (N != A.size(0) /*  || N != A.size(1) alow for rectangular matrices */)
    error("Incorrect dimension of matrix for application of boundary conditions."
          "Did you assemble it on a different mesh?");
  if (N != b.size())
    error("Incorrect dimension of matrix for application of boundary conditions."
        "Did you assemble it on a different mesh?");

  // A map to hold the mapping from boundary dofs to boundary values
  _map<uint, real> boundary_values;

  // Compute dofs and values
  computeBC(boundary_values, data);

  // Copy boundary value data to arrays
  uint* dofs = new uint[boundary_values.size()];
  real* values = new real[boundary_values.size()];
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
    real* x_values = new real[boundary_values.size()];
    x->get(x_values, boundary_values.size(), dofs);
    for (uint i = 0; i < boundary_values.size(); i++)
      values[i] -= x_values[i];
    delete[] x_values;
  }

  message("Applying boundary conditions to linear system.");

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
void DirichletBC::zero(GenericMatrix& A, Form const& form)
{
  // Get local data for application of boundary conditions
  BoundaryCondition::LocalData& data = this->updateLocalData(form);

  //
  DofMap const& dof_map = form.dofmaps()[1];

  // Simple check
  const uint N = dof_map.global_dimension();
  if (N != A.size(0))
    error("Incorrect dimension of matrix for application of boundary conditions."
          "Did you assemble it on a different mesh?");

  // A map to hold the mapping from boundary dofs to boundary values
  _map<uint, real> boundary_values;

  // Compute dofs and values
  computeBC(boundary_values, data);

  // Copy boundary value data to arrays
  uint* dofs = new uint[boundary_values.size()];
  _map<uint, real>::const_iterator boundary_value;
  uint i = 0;
  for (boundary_value = boundary_values.begin();
       boundary_value != boundary_values.end(); ++boundary_value)
    dofs[i++] = boundary_value->first;

  // Modify linear system (A_ii = 1)
  A.zero(boundary_values.size(), dofs);

  // Finalise changes to A
  A.apply();

  // Clear temporary arrays
  delete[] dofs;
}
//-----------------------------------------------------------------------------
void DirichletBC::initFromSubDomain(const SubDomain& sub_domain)
{
  dolfin_assert(facets_.size() == 0);

  // FIXME: This can be made more efficient, we should be able to
  // FIXME: extract the facets without first creating a MeshFunction on
  // FIXME: the entire mesh and then extracting the subset. This is done
  // FIXME: mainly for convenience (we may reuse mark() in SubDomain).

  // Make sure the mesh has been ordered
  mesh().order();

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
void DirichletBC::initFromMeshFunction(MeshFunction<uint>& sub_domains,
                                       uint sub_domain)
{
  dolfin_assert(facets_.size() == 0);

  // Make sure we have the facet - cell connectivity
  const uint dim = mesh().topology().dim();
  mesh().init(dim - 1, dim);

  // Make sure the mesh has been ordered
  mesh().order();

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
void DirichletBC::computeBC(_map<uint, real>& boundary_values,
                            BoundaryCondition::LocalData& data)
{
  // Choose strategy
    switch (method_)
    {
      case topological:
      computeBCTopological(boundary_values, data);
      break;
      case geometric:
      computeBCGeometric(boundary_values, data);
      break;
      case pointwise:
      computeBCPointwise(boundary_values, data);
      break;
      default:
      error("Unknown method for application of boundary conditions.");
      break;
    }
  }
//-----------------------------------------------------------------------------
void DirichletBC::computeBCTopological(_map<uint, real>& boundary_values,
                                       BoundaryCondition::LocalData& data)
{
  // Special case
    if (facets_.size() == 0 && dolfin::MPI::numProcesses() == 1)
    {
      warning("Found no facets matching domain for boundary condition.");
      return;
    }

    // Iterate over facets
#ifndef NO_PROGRESS_BAR
    Progress p("Computing Dirichlet boundary values, topological search", facets_.size());
#endif
    for (uint f = 0; f < facets_.size(); f++)
    {
      // Get cell number and local facet number
      uint cell_number = facets_[f].first;
      uint facet_number = facets_[f].second;

      // Create cell
      Cell cell(mesh(), cell_number);
      UFCCell ufc_cell(cell);

      ufc_cell.update(cell, mesh().distdata());

      // Interpolate function on cell
      g_.interpolate(data.w, ufc_cell, *data.finite_element, cell, facet_number);

      // Tabulate dofs on cell
      //data.dof_map->tabulate_dofs(data.cell_dofs, ufc_cell);
      data.dof_map->tabulate_dofs(data.cell_dofs, ufc_cell, cell.index());

      // Tabulate which dofs are on the facet
      data.dof_map->tabulate_facet_dofs(data.facet_dofs, facet_number);

      // Debugging print:
      /*
       cout << endl << "Handling BC's for:" << endl;
       cout << "Cell:  " << facet.entities(facet.dim() + 1)[0] << endl;
       cout << "Facet: " << local_facet << endl;
       */

      // Pick values for facet
      for (uint i = 0; i < data.dof_map->num_facet_dofs(); i++)
      {
        const uint dof = data.offset + data.cell_dofs[data.facet_dofs[i]];
        const real value = data.w[data.facet_dofs[i]];
        boundary_values[dof] = value;
        //cout << "Setting BC value: i = " << i << ", dof = " << dof << ", value = " << value << endl;
      }

#ifndef NO_PROGRESS_BAR
    p++;
#endif

  }
}
//-----------------------------------------------------------------------------
void DirichletBC::computeBCGeometric(_map<uint, real>& boundary_values,
                                     BoundaryCondition::LocalData& data)
{
  // Special case
    if (facets_.size() == 0 && dolfin::MPI::numProcesses() == 1)
    {
      warning("Found no facets matching domain for boundary condition.");
      return;
    }

    // Initialize facets, needed for geometric search
    message("Computing facets, needed for geometric application of boundary conditions.");
    uint const facet_dim = mesh().topology().dim() - 1;
    mesh().init(facet_dim);

    // Iterate over facets
#ifndef NO_PROGRESS_BAR
    Progress p("Computing Dirichlet boundary values, geometric search", facets_.size());
#endif
    CellType& celltype = mesh().type();
    Point dof_node;
    for (uint f = 0; f < facets_.size(); ++f)
    {
      // Get cell number and local facet number
      uint cell_number = facets_[f].first;
      uint facet_number = facets_[f].second;

      // Create facet
      Cell cell(mesh(), cell_number);
      Facet facet(mesh(), cell.entities(mesh().topology().dim() - 1)[facet_number]);

      // Loop the vertices associated with the facet
      for (VertexIterator vertex(facet); !vertex.end(); ++vertex)
      {
        // Loop the cells associated with the vertex
        for (CellIterator c(*vertex); !c.end(); ++c)
        {
          UFCCell ufc_cell(*c);
          ufc_cell.update(*c, mesh().distdata());
          bool interpolated = false;

          // Tabulate coordinates of dofs on cell
          data.dof_map->tabulate_coordinates(data.coordinates, ufc_cell);

          // Loop over all dofs on cell
          for (uint i = 0; i < data.dof_map->local_dimension(); ++i)
          {
            // Copy coordinates to node Point
            std::memcpy(&dof_node[0], &data.coordinates[i][0],
                        facet_dim*sizeof(real));
            // Check if the coordinates are on current facet and thus on boundary
            if (!celltype.intersects(facet, dof_node))
            {
              continue;
            }

            if(!interpolated)
            {

              // Tabulate dofs on cell
              data.dof_map->tabulate_dofs(data.cell_dofs, ufc_cell);

              // Interpolate function on cell
              g_.interpolate(data.w, ufc_cell, *data.finite_element, *c);
            }

            // Set boundary value
            const uint dof = data.offset + data.cell_dofs[i];
            const real value = data.w[i];
            boundary_values[dof] = value;
          }

        }
      }
    }
  }
//-----------------------------------------------------------------------------
void DirichletBC::computeBCPointwise(_map<uint, real>& boundary_values,
                                     BoundaryCondition::LocalData& data)
{
  // Iterate over cells
#ifndef NO_PROGRESS_BAR
    Progress p("Computing Dirichlet boundary values, pointwise search", mesh().numCells());
#endif
    for (CellIterator cell(mesh()); !cell.end(); ++cell)
    {
      UFCCell ufc_cell(*cell);
      ufc_cell.update(*cell, mesh().distdata());
      // Tabulate coordinates of dofs on cell
      data.dof_map->tabulate_coordinates(data.coordinates, ufc_cell);

      // Interpolate function only once and only on cells where necessary
      bool interpolated = false;

      // Loop all dofs on cell
      for (uint i = 0; i < data.dof_map->local_dimension(); ++i)
      {
        // Check if the coordinates are part of the sub domain
        if ( !this->sub_domain().inside(data.coordinates[i], false) )
          continue;

        if(!interpolated)
        {
          interpolated = true;
          // Tabulate dofs on cell
          data.dof_map->tabulate_dofs(data.cell_dofs, ufc_cell);
          // Interpolate function on cell
          g_.interpolate(data.w, ufc_cell, *data.finite_element, *cell);
        }

        // Set boundary value
        const uint dof = data.offset + data.cell_dofs[i];
        const real value = data.w[i];
        boundary_values[dof] = value;
      }

#ifndef NO_PROGRESS_BAR
      p++;
#endif
  }
}
//-----------------------------------------------------------------------------
void DirichletBC::getBC(uint n, uint* indicators, double* values,
                        DofMap const& dof_map, Form const& form)
{
  // Create local data for application of boundary conditions
  BoundaryCondition::LocalData& data = this->updateLocalData(form);

  // A map to hold the mapping from boundary dofs to boundary values
  _map<uint, real> boundary_values;

  // Compute dofs and values
  computeBC(boundary_values, data);

  if ( n != dof_map.global_dimension())
  {
    error("The n should be the same as dof_map.global_dimension()");
  }

  _map<uint, real>::const_iterator boundary_value;
  uint i = 0;
  for (boundary_value = boundary_values.begin();
      boundary_value != boundary_values.end(); ++boundary_value)
  {
    i = boundary_value->first;
    indicators[i] = 1;
    values[i] = boundary_value->second;
  }
}
//-----------------------------------------------------------------------------

}
