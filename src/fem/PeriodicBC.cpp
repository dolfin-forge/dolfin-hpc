// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007
//
// First added:  2007-07-08
// Last changed: 2008-04-22

#include <dolfin/fem/PeriodicBC.h>

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/UFCCell.h>

#include <vector>
#include <map>

using namespace dolfin;

// Comparison operator for hashing coordinates. Note that two
// coordinates are considered equal if equal to within round-off.
struct lt_coordinate
{
  bool operator()(const std::vector<real>& x, const std::vector<real>& y) const
  {
    if (x.size() > (y.size() + DOLFIN_EPS))
      return false;

    for (unsigned int i = 0; i < x.size(); i++)
    {
      if (x[i] < (y[i] - DOLFIN_EPS))
        return true;
      else if (x[i] > (y[i] + DOLFIN_EPS))
        return false;
    }

    return false;
  }
};

//-----------------------------------------------------------------------------
PeriodicBC::PeriodicBC(Mesh& mesh, const SubDomain& sub_domain) :
    BoundaryCondition("Periodic", mesh, sub_domain)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PeriodicBC::PeriodicBC(Mesh& mesh, const SubDomain& sub_domain,
                       const SubSystem& sub_system) :
    BoundaryCondition("Periodic", mesh, sub_domain, sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PeriodicBC::~PeriodicBC()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form)
{
  message("Applying periodic boundary conditions to linear system.");

  // FIXME: Make this work for non-scalar subsystems, like vector-valued
  // FIXME: Lagrange where more than one per element is associated with
  // FIXME: each coordinate. Note that globally there may very well be
  // FIXME: more than one dof per coordinate (for conforming elements).
  uint const gdim = mesh().geometry().dim();
  uint const tdim = mesh().geometry().dim();

  // Table of mappings from coordinates to dofs
  std::map<std::vector<real>, std::pair<int, int>, lt_coordinate> coordinate_dofs;
  typedef std::map<std::vector<real>, std::pair<int, int>, lt_coordinate>::iterator iterator;
  std::vector<real> xx(gdim);

  // Array used for mapping coordinates
  real * y = new real[gdim];
  for (uint i = 0; i < gdim; i++)
  {
    y[i] = 0.0;
  }

  // Make sure we have the facet - cell connectivity
  mesh().init(tdim - 1, tdim);

  // Create local data for application of boundary conditions
  FiniteElementSpace const& space = form.trial_space();
  ScratchSpace scratch(space);

  // Iterate over the facets of the mesh
#ifndef NO_PROGRESS_BAR
  Progress p("Applying periodic boundary conditions", mesh().size(tdim - 1));
#endif
  for (FacetIterator facet(mesh()); !facet.end(); ++facet)
  {
    // Get cell to which facet belongs (there may be two, but pick first)
    Cell cell(mesh(), facet->entities(tdim)[0]);
    scratch.cell.update(cell, mesh().distdata());

    // Get local index of facet with respect to the cell
    const uint local_facet = cell.index(*facet);

    // Tabulate dofs on cell
    scratch.dof_map->tabulate_dofs(scratch.dofs, scratch.mesh, scratch.cell);

    // Tabulate coordinates on cell
    scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);

    // Tabulate which dofs are on the facet
    scratch.dof_map->tabulate_facet_dofs(scratch.facet_dofs, local_facet);

    // Iterate over facet dofs
    for (uint i = 0; i < scratch.dof_map->num_facet_dofs(); ++i)
    {
      // Get dof and coordinate of dof
      uint const local_dof = scratch.facet_dofs[i];
      int const global_dof = static_cast<int>(scratch.offset + scratch.dofs[local_dof]);
      real *& x = scratch.coordinates[local_dof];

      // Map coordinate from H to G
      for (uint j = 0; j < gdim; j++)
      {
        y[j] = x[j];
      }
      sub_domain().map(x, y);

      // Check if coordinate is inside the domain G or in H
      const bool on_boundary = facet->numEntities(tdim) == 1;
      if (sub_domain().inside(x, on_boundary))
      {
        // Copy coordinate to std::vector
        for (uint j = 0; j < gdim; j++)
        {
          xx[j] = x[j];
        }

        // Check if coordinate exists from before
        iterator it = coordinate_dofs.find(xx);
        if (it != coordinate_dofs.end())
        {
          // Check that we don't have more than one dof per coordinate
          /*
           if (it->second.first != -1)
           {
           cout << "Coordinate: x =";
           for (uint j = 0; j < gdim; j++)
           cout << " " << xx[j];
           cout << endl;
           cout << "Degrees of freedom: " << it->second.first << " " << global_dof << endl;
           error("More than one dof associated with coordinate. Did you forget to specify the subsystem?");
           }
           */
          it->second.first = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal second value)
          std::pair<int, int> dofs(global_dof, -1);
          coordinate_dofs[xx] = dofs;
        }
      }
      else if (sub_domain().inside(y, on_boundary))
      {
        // y = F(x) is in G, so coordinate x is in H

        // Copy coordinate to std::vector
        for (uint j = 0; j < gdim; j++)
          xx[j] = y[j];

        // Check if coordinate exists from before
        iterator it = coordinate_dofs.find(xx);
        if (it != coordinate_dofs.end())
        {
          // Check that we don't have more than one dof per coordinate
          /*
           if (it->second.second != -1)
           {
           cout << "Coordinate: x =";
           for (uint j = 0; j < gdim; j++)
           cout << " " << xx[j];
           cout << endl;
           cout << "Degrees of freedom: " << it->second.second << " " << global_dof << endl;
           error("More than one dof associated with coordinate. Did you forget to specify the subsystem?");
           }
           */
          it->second.second = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal first value)
          std::pair<int, int> dofs(-1, global_dof);
          coordinate_dofs[xx] = dofs;
        }
      }
    }
#ifndef NO_PROGRESS_BAR
    p++;
#endif
  }

  // Insert 1 at (dof0, dof1)
  uint* rows = new uint[coordinate_dofs.size()];
  uint i = 0;
  for (iterator it = coordinate_dofs.begin(); it != coordinate_dofs.end(); ++it)
    rows[i++] = static_cast<uint>(it->second.first);
  A.ident(coordinate_dofs.size(), rows);

  // Insert -1 at (dof0, dof1) and 0 on right-hand side
  uint* cols = new uint[1];
  real* vals = new real[1];
  real* zero = new real[1];
  for (iterator it = coordinate_dofs.begin(); it != coordinate_dofs.end(); ++it)
  {
    // Check that we got both dofs
    const int dof0 = it->second.first;
    const int dof1 = it->second.second;
    if (dof0 == -1 || dof1 == -1)
    {
      cout << "At coordinate: x =";
      for (uint j = 0; j < gdim; j++)
        cout << " " << it->first[j];
      cout << endl;
      error(
          "Unable to find a pair of matching dofs for periodic boundary condition.");
    }

    // FIXME: Perhaps this can be done more efficiently?

    // Set x_i - x_j = 0
    rows[0] = static_cast<uint>(dof0);
    cols[0] = static_cast<uint>(dof1);
    vals[0] = -1;
    zero[0] = 0.0;

    A.set(vals, 1, rows, 1, cols);
    b.set(zero, 1, rows);
  }
  delete[] rows;
  delete[] cols;
  delete[] vals;
  delete[] zero;
  delete[] y;

  // Apply changes
  A.apply();
  b.apply();
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericMatrix& A, GenericVector& b,
                       const GenericVector& x, const BilinearForm& form)
{
  error("Periodic boundary conditions not implemented for nonlinear systems.");
}
//-----------------------------------------------------------------------------
