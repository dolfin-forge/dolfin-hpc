// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2015.
// Modified by Aurélien Larcher, 2013-2014. (rewrite, extension to any element)
//
// First added:  2007-05-01
// Last changed: 2014-05-22

#include <dolfin/fem/UFC.h>

#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SlipBC.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

#include <cmath>
#include <cstring>
#include <map>

#if (__sgi)
#define fmax(a,b) (a > b ? a : b) ;
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain const& sub_domain) :
    BoundaryCondition("SlipBC", mesh, sub_domain),
    mesh(mesh),
    node_normal(new NodeNormal(mesh, sub_domain)),
    node_normal_local(true),
    As(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain const& sub_domain, NodeNormal& normals) :
    BoundaryCondition("SlipBC", mesh, sub_domain),
    mesh(mesh),
    node_normal(&normals),
    node_normal_local(false),
    As(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain) :
    BoundaryCondition("SlipBC", sub_domains, sub_domain),
    mesh(sub_domains.mesh()),
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain const& sub_domain,
               SubSystem const& sub_system) :
    BoundaryCondition("SlipBC", mesh, sub_domain, sub_system),
    mesh(mesh),
    node_normal(new NodeNormal(mesh, sub_domain)),
    node_normal_local(true),
    As(NULL)
{
  // Do nothing
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
               SubSystem const& sub_system) :
    BoundaryCondition("SlipBC", sub_domains, sub_domain, sub_system),
    mesh(sub_domains.mesh()),
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::~SlipBC()
{
  if (node_normal_local) delete node_normal;

  //
  delete As;
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form)
{
  FiniteElementSpace const& fullspace = form.test_space();

  // If the subsystem is not empty then the scratch space corresponds to the
  // given subspace, otherwise it is the space itself.
  ScratchSpace scratch(fullspace, sub_system());
  // Optimize loops if the space is linear Lagrange: the space dimension is
  // the number of vertices times the number of components
  bool const is_P1 = (scratch.space_dimension
                        == mesh.type().numEntities(0) * scratch.size);
  
  const std::string la_backend = dolfin_get("linear algebra backend");
    
  if (As == NULL || As->size(0) != A.size(0) || As->size(1) != A.size(1))
  {
    // Create data structure for local assembly data
    if (la_backend == "JANPACK")
    {
      As = static_cast<Matrix *>(&A);
    }
    else
    {
      delete As;
      As = new Matrix();
#if HAVE_PETSC
      (*(As->instance())).down_cast<PETScMatrix>().dup(A);
#endif
    }

    // Initialize ghosts for rhs vector using the full space dofmap
    std::set<uint> rows;
    std::map<uint, uint> mapping;
    DofMap const& fulldofmap = fullspace.dofmap();
    uint * celldofs = new uint[fulldofmap.local_dimension()];
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      scratch.cell.update(*c, mesh.distdata());
      fulldofmap.tabulate_dofs(celldofs, scratch.cell, *c);
      for (uint j = 0; j < fulldofmap.local_dimension(); ++j)
      {
        rows.insert(celldofs[j]);
      }
    }
    delete[] celldofs;
    b.init_ghosted(rows.size(), rows, mapping);

    // Initialize normal field on given space and compute at the boundary
    if (sub_system().depth() == 0)
    {
      node_normal->init(fullspace);
    }
    else
    {
      FiniteElementSpace subspace(fullspace, sub_system());
      node_normal->init(subspace);
    }
    node_normal->compute();

    // Create boundary markers for given topological dimension if the subdomain
    // is defined geometrically.
    if (this->has_geometrical_sub_domain())
    {
      if (is_P1)
      {
        // Markers are vertex-based
        BoundaryCondition::init_markers(0);
      }
      else
      {
        // Markers are facet based
        BoundaryCondition::init_markers(mesh.topology().dim() - 1);
      }
    }

    // Initialize local data structures
    std::set<uint>::iterator it = row_indices.begin();
    for (uint i = 0; i < mesh.type().dim(); ++i)
    {
      a[i].reserve(A.size(0));
      a_slip_row[i].reserve(A.size(0));
      a_col_indices[i].reserve(A.size(0));
      it = row_indices.insert(it, i);
    }
  }

  // Copy global stiffness matrix into temporary one (if needed)
  if (la_backend != "JANPACK")
  {
    *(As->instance()) = A;
  }


  // Use legacy vertex-based implementation if Lagrange P1.
  if (is_P1)
  {
    applySlipBC_P1(A, b, form, scratch);
  }
  else
  {
    applySlipBC(A, b, form, scratch);
  }

  // Apply changes in the temporary matrix
  As->apply();

  // Apply changes in the stiffness matrix and load vector (if needed)
  if (la_backend != "JANPACK")
  {
    A = *(As->instance());
  }
  b.apply();

}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
                   BilinearForm const& form)
{
  error("SlipBC not implemented for non linear systems");
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC_P1(GenericMatrix& A, GenericVector& b,
                            BilinearForm const& form, ScratchSpace& scratch)
{
  BoundaryMesh& boundary = mesh.exterior_boundary();
  if (boundary.numCells())
  {
    MeshFunction<uint> const& sub_domains = this->sub_domain_markers();
    uint sub_domain_idx = this->sub_domain_index();

    // Used to get the global dof indices
    DofMap const& Udofmap = form.test_space().dofmap();
    uint * fulldofs = new uint[Udofmap.local_dimension()];

    // If any, to get the node id
    bool const same_space = (sub_system().depth() == 0);
    DofMap const& Ndofmap = node_normal->basis()[0].space().dofmap();

    Array<uint> node_Udofs;
    Array<uint> node_Ndofs;
    uint const gdim = mesh.type().dim();
    uint const vdim = mesh.type().numEntities(0); // number of nodes = vertices

    for (VertexIterator v(boundary); !v.end(); ++v)
    {
      Vertex vertex(mesh, boundary.vertex_index(*v));

      // Skip vertices not inside the sub domain
      if (sub_domains(vertex) != sub_domain_idx)
      {
        continue;
      }

      uint node = vertex.index();
      if (!mesh.distdata().is_ghost(node, 0))
      {
        Cell cell(mesh, (vertex.entities(gdim))[0]);
        scratch.cell.update(cell, mesh.distdata());

        // Find the vertex position in the cell
        uint *cvi = cell.entities(0);
        uint ci = 0;
        for (ci = 0; ci < cell.numEntities(0); ++ci)
        {
          if (cvi[ci] == node) break;
        }
        uint const vindex = ci;

        // Get the component dofs of U for the given node
        uint voff = 0;
        Udofmap.tabulate_dofs(fulldofs, scratch.cell, cell);
        for (uint i = 0; i < scratch.size; ++i, voff += vdim)
        {
          node_Udofs.push_back(fulldofs[scratch.offset + voff + vindex]);
        }

        // The node id is the global index of the 1st comp. of the node normal
        // which is the same as the fullspace dof is there is no subsystem
        if (same_space)
        {
          uint node_id = node_Udofs[0];
          applyNodeBC(A, b, mesh, node_id, node_Udofs, node_Udofs);
          node_Udofs.clear();
        }
        else
        {
          // Get the component dofs of N for the given node
          uint voff = 0;
          Ndofmap.tabulate_dofs(scratch.dofs, scratch.cell, cell);
          for (uint i = 0; i < scratch.size; ++i, voff += vdim)
          {
            node_Ndofs.push_back(scratch.dofs[voff + vindex]);
          }
          uint node_id = node_Ndofs[0];
          applyNodeBC(A, b, mesh, node_id, node_Udofs, node_Ndofs);
          node_Udofs.clear();
          node_Ndofs.clear();
        }
      }
    }

    delete[] fulldofs;
  }
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC(GenericMatrix& A, GenericVector& b,
                         BilinearForm const& form, ScratchSpace& scratch)
{
  // Be careful for now
  uint const gdim = mesh.geometry().dim();
  dolfin_assert(scratch.size == gdim);

  BoundaryMesh& boundary = mesh.exterior_boundary();
  if (boundary.numCells())
  {
    MeshFunction<uint> const& sub_domains = this->sub_domain_markers();
    uint const sub_domain_idx = this->sub_domain_index();

    DofMap const& Udofmap = form.test_space().dofmap();
    uint * fulldofs = new uint[Udofmap.local_dimension()];

    // If any, to get the node id defined as the dof of the first normal comp.
    bool const same_space = (sub_system().depth() == 0);
    DofMap const& Ndofmap = node_normal->basis()[0].space().dofmap();

    //
    Array<uint> node_Udofs;
    Array<uint> node_Ndofs;
    uint const num_facet_dofs = scratch.dof_map->num_facet_dofs();
    uint const num_facet_nodes = num_facet_dofs / scratch.size;
    uint const vdim = scratch.local_dimension / scratch.size; // component offset
    _set<uint> visited_nodes;
    for (CellIterator c(boundary); !c.end(); ++c)
    {
      Facet facet(mesh, boundary.facet_index(*c));

      // Skip facets outside the sub domain
      if (sub_domains(facet) != sub_domain_idx)
      {
        continue;
      }

      //
      Cell cell(mesh, facet.entities(gdim)[0]);
      uint const local_facet = cell.index(facet);
      scratch.cell.update(cell, mesh.distdata());
      scratch.dof_map->tabulate_coordinates(scratch.coordinates, scratch.cell);
      // Attention, local-to-local mapping.
      scratch.dof_map->tabulate_facet_dofs(scratch.facet_dofs, local_facet);

      // Tabulate full space and nodenormal space if needed
      Udofmap.tabulate_dofs(fulldofs, scratch.cell, cell);
      if (!same_space)
      {
        Ndofmap.tabulate_dofs(scratch.dofs, scratch.cell, cell);
      }

      // Apply slipbc on each non ghosted facet dof node of the subspace
      for (uint ni = 0; ni < num_facet_nodes; ++ni)
      {
        // Get the cell local index of the first dof of node ni
        uint const ni_celldof0 = scratch.facet_dofs[ni];
        // Get the global index of the first dof of node ni
        uint const ni_Udof0 = fulldofs[scratch.offset + ni_celldof0];

        // Skip the node if it has already been encountered
        if (visited_nodes.count(ni_Udof0))
        {
          continue;
        }
        else
        {
          visited_nodes.insert(ni_Udof0); // mark as visited
        }
        // Skip the node if its component dofs are ghosted
        if (Udofmap.is_ghost(ni_Udof0))
        {
          continue;
        }
        // Skip the node if the subdomain is defined geometrically
        if (this->has_geometrical_sub_domain()
            && !sub_domain().inside(scratch.coordinates[ni_celldof0], true))
        {
          continue;
        }

        // Get the component dofs of U for the given node
        uint voff = 0;
        for (uint i = 0; i < scratch.size; ++i, voff += vdim)
        {
          uint ii = fulldofs[scratch.offset + voff + ni_celldof0];
          node_Udofs.push_back(ii);
          dolfin_assert(!Udofmap.is_ghost(ii));
        }

        // The node id is the global index of the 1st comp. of the node normal
        // which is the same as the fullspace dof is there is no subsystem
        if (same_space)
        {
          uint node_id = node_Udofs[0];
          applyNodeBC(A, b, mesh, node_id, node_Udofs, node_Udofs);
          node_Udofs.clear();
        }
        else
        {
          // Get the component dofs of N for the given node
          uint voff = 0;
          for (uint i = 0; i < scratch.size; ++i, voff += vdim)
          {
            uint ii = scratch.dofs[voff + ni_celldof0];
            node_Ndofs.push_back(ii);
          }
          uint node_id = node_Ndofs[0];
          applyNodeBC(A, b, mesh, node_id, node_Udofs, node_Ndofs);
          node_Udofs.clear();
          node_Ndofs.clear();
        }
      }
    }
    visited_nodes.clear();
    delete[] fulldofs;
  }
}

//-----------------------------------------------------------------------------
void SlipBC::applyNodeBC(GenericMatrix& A, GenericVector& b, Mesh const& mesh,
                         uint const node_id, Array<uint> const& Udofs,
                         Array<uint> const& Ndofs)
{
  // Naive reimplementation -- Aurélien
  // The node type defines the number of discriminated surfaces at the node.
  // Therefore it is the number of constrained directions up to the topological
  // dimension
  uint const gdim = mesh.topology().dim();
  uint const n_type = std::min(node_normal->node_type(node_id), gdim);
  dolfin_assert(n_type > 0);

  // Initialize set of row indices for reordering
  std::set<uint> row_idx = row_indices;
  std::set<uint>::iterator it = row_idx.begin();

  //--- Fill data structures for each space coordinate ---
  Array<Function>& basis_functions = node_normal->basis();
  for (uint i = 0; i < gdim; ++i)
  {
    // Copy non-zero entries from the stiffness matrix into local row
    A.getrow(Udofs[i], a_col_indices[i], a[i]);

    // Reset rhs for slip
    uint nb_cols = a_col_indices[i].size();
    a_slip_row[i].resize(nb_cols);
    std::fill(a_slip_row[i].begin(), a_slip_row[i].end(), 0.0);

    // Reset lhs for slip
    l_slip[i] = 0.0;

    // Zero the row in the matrix copy
    As->set(&a_slip_row[i][0], 1, &Udofs[i], nb_cols, &a_col_indices[i][0]);

    // Fill component i-th basis vector
    real (&v)[EuclideanSpace::MAX_DIMENSION] = basis_[i];
    basis_functions[i].vector().get(&v[0], gdim, &Ndofs[0]);

    // Determine maximum component (to be simplified)
    it = row_idx.begin();
    max[i] = (*it);
    for (; it != row_idx.end(); ++it)
    {
      if (std::fabs(v[*it]) > std::fabs(v[max[i]]))
      {
        max[i] = *it;
      }
    }
    row_idx.erase(max[i]);
  }

  // Integer n_type represents the number of constrained directions
  if (n_type < gdim)  // At least one free direction
  {
    //--- For each constrained direction
    for (uint i = 0; i < n_type; ++i)
    {
      // Set row to the global index
      row[i] = Udofs[max[i]];

      // Update the LHS row with the vector components
      As->set(&basis_[i][0], 1, &row[i], gdim, &Udofs[0]);

      // Reset rhs for slip
      l_slip[i] = 0.0;
    }
    //--- Copy RHS to local vector (here since a surface node is most probable)
    b.get(&l[0], gdim, &Udofs[0]);
    //--- TODO: Row swapping of free directions assumes same row structure
    if (gdim == 3 && n_type == 1)
    {
      // Find position of the diagonal columns in the vectors
      uint const maxn_1 = max[1];
      uint diag_idx_1 = 0;
      while (a_col_indices[maxn_1][diag_idx_1] != Udofs[maxn_1])
      {
        ++diag_idx_1;
      }
      real row_d1 = std::fabs(
          a[0][diag_idx_1] * basis_[1][0] + a[1][diag_idx_1] * basis_[1][1]
              + a[2][diag_idx_1] * basis_[1][2])
          - std::fabs(
              a[0][diag_idx_1] * basis_[2][0] + a[1][diag_idx_1] * basis_[2][1]
                  + a[2][diag_idx_1] * basis_[2][2]);

      uint const maxn_2 = max[2];
      uint diag_idx_2 = 0;
      while (a_col_indices[maxn_2][diag_idx_2] != Udofs[maxn_2])
      {
        ++diag_idx_2;
      }
      real row_d2 = std::fabs(
          a[0][diag_idx_2] * basis_[2][0] + a[1][diag_idx_2] * basis_[2][1]
              + a[2][diag_idx_2] * basis_[2][2])
          - std::fabs(
              a[0][diag_idx_2] * basis_[1][0] + a[1][diag_idx_2] * basis_[1][1]
                  + a[2][diag_idx_2] * basis_[1][2]);

      if (row_d1 < 0.0 || row_d2 < 0.0)
      {
        max[1] = maxn_2;
        max[2] = maxn_1;
      }
    }
    //--- For each free direction
    for (uint i = n_type; i < gdim; ++i)
    {
      // Set row to the global index
      row[i] = Udofs[max[i]];

      // Project equation on the tangential vector: a[j].tau_i
      // Beware, we assume here that the component contributions were
      // inserted in the same order !
      uint const nb_cols = a_col_indices[max[i]].size();
      for (uint k = 0; k < gdim; ++k)
      {
        for (uint j = 0; j < nb_cols; ++j)
        {
          a_slip_row[i][j] += a[k][j] * basis_[i][k];
        }
        l_slip[i] += l[k] * basis_[i][k];
      }
      // Update the LHS row with the projection of the equation on tau_i
      As->setrow(row[i], a_col_indices[max[i]], a_slip_row[i]);
    }
    //--- Apply local equation RHS to the copy of the matrix
    b.set(&l_slip[0], gdim, &row[0]);
  }
  else // All directions are constrained
  {
    for (uint i = 0; i < gdim; ++i)
    {
      // Find position of the diagonal in the vectors
      uint diag_idx = 0;
      while (a_col_indices[i][diag_idx] != Udofs[i])
      {
        ++diag_idx;
      }

      // Scale the diagonal entry and update the LHS diagonal [Disabled]
      real diag_val = 1.; //std::fabs(a[i][diag_idx]);
      As->set(&diag_val, 1, &Udofs[i], 1, &Udofs[i]);
    }
    //Apply local equation RHS to the copy of the matrix
    b.set(&l_slip[0], gdim, &Udofs[0]);
  }

  //--- DEBUGGING
  //real entries[3] = { 1.0, 2.0, 3.0 };
  //b.set(&entries[0], gdim, &Udofs[0]);
  //
}
//-----------------------------------------------------------------------------

}

