// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2012.
//
// First added:  2007-05-01
// Last changed: 2012-03-04

#include <dolfin/fem/UFC.h>

#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/fem/SlipBC.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/main/MPI.h>
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
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(0)
{
  // Initialize sub domain markers on vertices
  BoundaryCondition::init_markers(0);
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain const& sub_domain, NodeNormal& normals) :
    BoundaryCondition("SlipBC", mesh, sub_domain),
    mesh(mesh),
    node_normal(&normals),
    node_normal_local(false),
    As(0)
{
  // Initialize sub domain markers
  BoundaryCondition::init_markers(0);
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain) :
    BoundaryCondition("SlipBC", sub_domains, sub_domain),
    mesh(sub_domains.mesh()),
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain const& sub_domain,
               const SubSystem& sub_system) :
    BoundaryCondition("SlipBC", mesh, sub_domain, sub_system),
    mesh(mesh),
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(0)
{
  // Set sub domain markers
  BoundaryCondition::init_markers(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
               const SubSystem& sub_system) :
    BoundaryCondition("SlipBC", sub_domains, sub_domain, sub_system),
    mesh(sub_domains.mesh()),
    node_normal(new NodeNormal(mesh)),
    node_normal_local(true),
    As(0)
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
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form)
{
  FiniteElementSpace const& space = form.test_space();
  ScratchSpace scratch(space);

  if (As == 0)
  {
    // Create data structure for local assembly data
    const std::string la_backend = dolfin_get("linear algebra backend");
    if (la_backend == "JANPACK")
    {
      As = new Matrix(A.size(0), A.size(1));
      *(As->instance()) = A;
      //      (*(As->instance())).down_cast<JANPACKMat>().dup(A);
    }
    else
    {
      As = new Matrix();
      (*(As->instance())).down_cast<PETScMatrix>().dup(A);
    }

    if (MPI::numProcesses() > 1)
    {
      std::map<uint, uint> mapping;
      for (CellIterator c(mesh); !c.end(); ++c)
      {
        scratch.cell.update(*c, mesh.distdata());
        space.dofmap().tabulate_dofs(scratch.dofs, scratch.cell, c->index());

        for (uint j = 0; j < scratch.local_dimension; ++j)
        {
          off_proc_rows.insert(scratch.dofs[j]);
        }
      }
      b.init_ghosted(off_proc_rows.size(), off_proc_rows, mapping);
    }

    // Initialize and compute normal field at the boundary
    node_normal->init(mesh, space.element().signature());
    node_normal->compute();

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

  // Copy global stiffness matrix into temporary one
  *(As->instance()) = A;

  // Use vertex-based implementation if Lagrange P1.
  if (space.element().space_dimension()
      == mesh.type().numEntities(0) * mesh.type().dim())
  {
    applySlipBC_P1(A, b, form, scratch);
  }
  else
  {
    error("SlipBC is not supported for other spaces than Lagrange P1");
  }

  // Apply changes in the temporary matrix
  As->apply();

  // Apply changes in the stiffness matrix and load vector
  A = *(As->instance());
  b.apply();

}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                   const BilinearForm& form)
{
  error("SlipBC not implemented for non linear systems");
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC_P1(GenericMatrix& A, GenericVector& b,
                            const BilinearForm& form, ScratchSpace& scratch)
{
  BoundaryMesh& boundary = mesh.exterior_boundary();
  MeshFunction<uint> * vertex_map = boundary.data().meshFunction("vertex map");
  if (boundary.numCells())
  {
    MeshFunction<uint> const& sub_domains = this->sub_domain_markers();
    uint sub_domain_idx = this->sub_domain_index();

    DofMap const& dofmap = form.test_space().dofmap();

    Array<uint> nodes;
    uint gdim = mesh.geometry().dim();
    uint cdim = mesh.type().numVertices(mesh.topology().dim());

    for (VertexIterator v(boundary); !v.end(); ++v)
    {

      Vertex vertex(mesh, vertex_map->get(*v));

      // Skip vertices not inside the sub domain
      if (sub_domains(vertex) != sub_domain_idx) continue;

      uint node = vertex.index();
      if (!mesh.distdata().is_ghost(node, 0) || MPI::numProcesses() == 1)
      {
        Cell cell(mesh, (vertex.entities(gdim))[0]);

        uint *cvi = cell.entities(0);
        uint ci = 0;
        for (ci = 0; ci < cell.numEntities(0); ++ci)
        {
          if (cvi[ci] == node) break;
        }

        scratch.cell.update(cell, mesh.distdata());
        dofmap.tabulate_dofs(scratch.dofs, scratch.cell, cell.index());

        for (uint i = 0; i < gdim; i++, ci += cdim)
        {
          nodes.push_back(scratch.dofs[ci]);
        }

        applyNodeBC(A, b, mesh, node, nodes);
        nodes.clear();
      }
    }
  }
}

//-----------------------------------------------------------------------------
void SlipBC::applyNodeBC(GenericMatrix& A, GenericVector& b, Mesh const& mesh,
                       uint const node, Array<uint> const& dofs)
{
  // Naive reimplementation -- Aurélien
  // The node type defines the number of discriminated surfaces at the node.
  // Therefore it is the number of constrained directions up to the topological
  // dimension
  real node_type = 0;
  node_type = node_normal->vertex_type.get(node);
  //node_normal->node_type().vector().get(&node_type, 1, &node);
  uint const tdim = mesh.topology().dim();
  int const n_type = std::min((int) std::floor(node_type), (int) tdim);

  // Initialize set of row indices for reordering
  std::set<uint> row_idx = row_indices;
  std::set<uint>::iterator it = row_idx.begin();

  //--- Fill data structures for each space coordinate ---
  Array<Function>& basis_functions = node_normal->basis();
  for (uint i = 0; i < tdim; ++i)
  {
    // Copy non-zero entries from the stiffness matrix into local row
    A.getrow(dofs[i], a_col_indices[i], a[i]);

    // Reset rhs for slip
    uint nb_cols = a_col_indices[i].size();
    a_slip_row[i].resize(nb_cols);
    std::fill(a_slip_row[i].begin(), a_slip_row[i].end(), 0.0);

    // Reset lhs for slip
    l_slip[i] = 0.0;

    // Zero the row in the matrix copy
    As->set(&a_slip_row[i][0], 1, &dofs[i], nb_cols, &a_col_indices[i][0]);

    // Fill component i-th basis vector
    real (&v)[3] = basis_[i];
    basis_functions[i].vector().get(&v[0], tdim, &dofs[0]);

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
  if (n_type < (int) tdim)  // At least one free direction
  {
    //--- For each constrained direction
    for (uint i = 0; i < n_type; ++i)
    {
      // Set row to the global index
      row[i] = dofs[max[i]];

      // Update the LHS row with the vector components
      As->set(&basis_[i][0], 1, &row[i], tdim, &dofs[0]);

      // Reset rhs for slip
      l_slip[i] = 0.0;
    }
    //--- Copy RHS to local vector (here since a surface node is most probable)
    b.get(&l[0], tdim, &dofs[0]);
    //--- For each free direction
    for (uint i = n_type; i < tdim; ++i)
    {
      // Set row to the global index
      row[i] = dofs[max[i]];

      // Project equation on the tangential vector: a[j].tau_i
      // Beware, we assume here that the component contributions were
      // inserted in the same order !
      uint const nb_cols = a_col_indices[max[i]].size();
      for (uint k = 0; k < tdim; ++k)
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
    b.set(&l_slip[0], tdim, &row[0]);
  }
  else // All directions are constrained
  {
    for (uint i = 0; i < tdim; ++i)
    {
      // Find position of the diagonal in the vectors
      uint diag_idx = 0;
      while (a_col_indices[i][diag_idx] != dofs[i])
      {
        ++diag_idx;
      }

      // Scale the diagonal entry and update the LHS diagonal
      real diag_val = std::fabs(a[i][diag_idx]);
      As->set(&diag_val, 1, &dofs[i], 1, &dofs[i]);
    }
    //Apply local equation RHS to the copy of the matrix
    b.set(&l_slip[0], tdim, &dofs[0]);
  }
}
//-----------------------------------------------------------------------------

}

