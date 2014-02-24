// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2012.
//
// First added:  2007-05-01
// Last changed: 2012-03-04

#include <dolfin/fem/SlipBC.h>

#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexNormal.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/parameter/parameters.h>

#include <ufc.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

#if (__sgi)
#define fmax(a,b) (a > b ? a : b) ;
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, const SubDomain& sub_domain) :
    BoundaryCondition("Slip"),
    mesh_(mesh),
    boundary_(mesh.exterior_boundary()),
    sub_domains_(NULL),
    sub_domain_(0),
    local_sub_domains_(false),
    user_sub_domain_(&sub_domain),
    normal_(new NodeNormal(mesh_, VertexNormal::none)),
    As_(NULL),
    tdim_(mesh.topology().dim())
{
  // Initialize sub domain markers
  init(sub_domain);

  sub_system_ = SubSystem(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(BoundaryNormal& normal, const SubDomain& sub_domain) :
    BoundaryCondition("Slip"),
    mesh_(normal.mesh()),
    boundary_(normal.boundary()),
    sub_domains_(NULL),
    sub_domain_(0),
    local_sub_domains_(false),
    user_sub_domain_(&sub_domain),
    local_normal_(false),
    normal_(&normal),
    As_(NULL),
    tdim_(mesh_.topology().dim())
{
  // Initialize sub domain markers
  init(sub_domain);

  sub_system_ = SubSystem(0);
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain) :
    BoundaryCondition("Slip"),
    mesh_(sub_domains.mesh()),
    boundary_(sub_domains.mesh().exterior_boundary()),
    sub_domains_(&sub_domains),
    sub_domain_(sub_domain),
    local_sub_domains_(false),
    sub_system_(),
    user_sub_domain_(NULL),
    local_normal_(true),
    normal_(new NodeNormal(mesh_, VertexNormal::none)),
    As_(NULL),
    N_local_(0),
    N_offset_(0),
    tdim_(mesh_.topology().dim())
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, const SubDomain& sub_domain,
               const SubSystem& sub_system) :
    BoundaryCondition("Slip"),
    mesh_(mesh),
    boundary_(mesh.exterior_boundary()),
    sub_domains_(NULL),
    sub_domain_(0),
    local_sub_domains_(false),
    sub_system_(sub_system),
    user_sub_domain_(&sub_domain),
    local_normal_(true),
    normal_(new NodeNormal(mesh_, VertexNormal::none)),
    As_(NULL),
    N_local_(0),
    N_offset_(0),
    tdim_(mesh.topology().dim())
{
  // Set sub domain markers
  init(sub_domain);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
               const SubSystem& sub_system) :
    BoundaryCondition("Slip"),
    mesh_(sub_domains.mesh()),
    boundary_(sub_domains.mesh().exterior_boundary()),
    sub_domains_(&sub_domains),
    sub_domain_(sub_domain),
    local_sub_domains_(false),
    sub_system_(sub_system),
    user_sub_domain_(NULL),
    local_normal_(true),
    normal_(new NodeNormal(mesh_, VertexNormal::none)),
    As_(NULL),
    N_local_(0),
    N_offset_(0),
    tdim_(mesh_.topology().dim())
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::~SlipBC()
{
  if (local_normal_)
  {
    delete normal_;
  }
  // Delete sub domain markers if created locally
  if (local_sub_domains_)
  {
    delete sub_domains_;
  }
  delete As_;
}
//-----------------------------------------------------------------------------
BoundaryNormal& SlipBC::normal()
{
  return *normal_;
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const Form& form)
{
  apply(A, b, form.dofmaps()[1], form);
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map,
                   const ufc::form& ufc_form)
{
  dolfin::error("Not implemented:\n",
                "void apply(GenericMatrix& A, GenericVector& b,\n",
                "\tDofMap const& dof_map, const ufc::form& form)");
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                   const Form& form)
{
  apply(A, b, form.dofmaps()[1], form);
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                   const DofMap& dof_map, const ufc::form& ufc_form)
{
  dolfin::error(
      "Not implemented:\n",
      "void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,\n",
      "\tDofMap const& dof_map, const ufc::form& form)");
}

//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map,
                   const Form& form)
{

  if (MPI::processNumber() == 0)
  {
    dolfin_set("output destination", "terminal");
    message("Applying SlipBC boundary conditions to linear system.");
    dolfin_set("output destination", "silent");
  }

  UFC ufc(form.form(), mesh_, form.dofmaps());

  if (As_ == NULL)
  {
    // Need to recompute the boundary normals.

    // Create data structure for local assembly data
    std::string const la_backend = dolfin_get("linear algebra backend");
    if (la_backend == "JANPACK")
    {
      As_ = new Matrix(A.size(0), A.size(1));
      *(As_->instance()) = A;
      //      (*(As->instance())).down_cast<JANPACKMat>().dup(A);
    }
    else
    {
      As_ = new Matrix();
      (*(As_->instance())).down_cast<PETScMatrix>().dup(A);
    }

    if (MPI::numProcesses() > 1)
    {
      std::map<uint, uint> mapping;
      for (CellIterator c(mesh_); !c.end(); ++c)
      {
        ufc.update(*c, mesh_.distdata());
        (form.dofmaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, c->index());

        for (uint j = 0; j < (form.dofmaps())[0].local_dimension(); ++j)
        {
          off_proc_rows_.insert(ufc.dofs[0][j]);
        }
      }

      b.init_ghosted(off_proc_rows_.size(), off_proc_rows_, mapping);
    }
  }
  // Copy global stiffness matrix into temporary one
  *(As_->instance()) = A;

  Array<uint> nodes;
  uint gdim = mesh_.geometry().dim();
  uint cdim = mesh_.type().numVertices(mesh_.topology().dim());

  uint count = 0;
  if (boundary_.numCells())
  {
    MeshFunction<uint> * vertex_map = boundary_.data().meshFunction(
        "vertex map");
    for (VertexIterator v(boundary_); !v.end(); ++v)
    {
      Vertex vertex(mesh_, vertex_map->get(*v));

      // Skip facets not inside the sub domain
      if ((*sub_domains_)(vertex) != sub_domain_)
      {
        continue;
      }

      uint node = vertex.index();
      if (!mesh_.distdata().is_ghost(node, 0) || MPI::numProcesses() == 1)
      {
        Cell cell(mesh_, (vertex.entities(gdim))[0]);

        uint *cvi = cell.entities(0);
        uint ci = 0;
        for (ci = 0; ci < cell.numEntities(0); ci++)
          if (cvi[ci] == node) break;

        ufc.update(cell, mesh_.distdata());
        (form.dofmaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, cell.index());

        // Get components of the vector-valued function at the current node.
        for (uint i = 0; i < gdim; i++, ci += cdim)
          nodes.push_back(ufc.dofs[0][ci]);

        applySlipBC((Matrix&) A, *As_, (Vector&) b, mesh_, node, nodes);
        count++;
        nodes.clear();
      }
    }
  }

  // Apply changes in the temporary matrix
  As_->apply();

  // Apply changes in the stiffness matrix and load vector
  A = *(As_->instance());
  b.apply();

}
//-----------------------------------------------------------------------------
void SlipBC::init(const SubDomain& sub_domain)
{
  // Create mesh function for sub domain markers on facets
  mesh_.init(0);
  sub_domains_ = new MeshFunction<uint>(mesh_, 0);
  local_sub_domains_ = true;

  // Mark everything as sub domain 1
  (*sub_domains_) = 1;

  // Mark the sub domain as sub domain 0
  sub_domain.mark(*sub_domains_, 0);
}

//-----------------------------------------------------------------------------
void SlipBC::applySlipBC(Matrix& A, Matrix& As, Vector& b, Mesh const& mesh,
                         uint const& node, Array<uint> const& dofs)
{
  // Naive reimplementation -- Aurélien

  // The node type defines the number of discriminated surfaces at the node.
  // Therefore it is the number of constrained directions up to the topological
  // dimension
  real node_type = 0;
  normal_->node_type().vector().get(&node_type, 1, &node);
  int const n_type = std::max((int) std::floor(node_type), (int) tdim_);

  // Initialize set of row indices for reordering
  std::set<uint> row_idx;
  std::set<uint>::iterator it = row_idx.begin();
  for (uint i = 0; i < tdim_; ++i)
  {
    it = row_idx.insert(it, i);
  }

  //--- Fill data structures for each space coordinate ---
  Array<Function>& basis_functions = normal_->basis();
  for (uint i = 0; i < tdim_; ++i)
  {
    // Copy non-zero entries from the stiffness matrix into local row
    A.getrow(dofs[i], a_col_indices[i], a[i]);

    // Copy rhs to local vector
    l[i] = b[dofs[i]];

    // Fill component i-th basis vector
    real (&v)[3] = basis_[i];
    basis_functions[i].vector().get(&v[0], tdim_, &dofs[0]);

    // Determine maximum component:
    for (it = row_idx.begin(); it != row_idx.end(); ++it)
    {
      if (std::fabs(v[*it]) > std::fabs(v[row[i]]))
      {
        row[i] = *it;
      }
    }
    row_idx.erase(row[i]);
  }

  //--- Apply local equation LHS to the copy of the matrix ---
  // n_type represents the number of constrained directions
  if (n_type < (int) tdim_)  // At least one free direction
  {

    // For each constrained direction
    for (uint i = 0; i < n_type; ++i)
    {
      // Set row to the global index
      row[i] = dofs[i];

      // Zero the row
      uint nb_cols = a_col_indices[i].size();
      a_slip_row.resize(nb_cols, 0.0);
      As.set(&a_slip_row[0], 1, &row[i], nb_cols, &a_col_indices[i][0]);

      As.zero(1, &row[i]);

      // Update the LHS row with the vector components
      As.set(&basis_[i][0], 1, &row[i], tdim_, &dofs[i]);

      // Reset rhs for slip
      l_slip[i] = 0.0;
    }

    // For each free direction
    for (uint i = n_type; i < tdim_; ++i)
    {
      // Set row to the global index
      row[i] = dofs[i];

      // Initialize row for application of the boundary condition
      uint nb_cols = a_col_indices[i].size();
      a_slip_row.resize(nb_cols, 0.0);

      // Project equation on the tangential vector: a[j].tau_i
      // Beware, we assume here that the component contributions were
      // inserted in the same order !
      l_slip[i] = 0.0;
      for (uint d = 0; d < tdim_; ++d)
      {
        for (uint j = 0; j < nb_cols; ++j)
        {
          a_slip_row[j] += a[d][j] * basis_[i][d];
        }
        l_slip[i] += l[d] * basis_[i][d];
      }

      // Update the LHS row with the projection of the equation on tau_i
      As.setrow(row[i], a_col_indices[i], a_slip_row);
    }

  }
  else // All directions are constrained
  {
    for (uint i = 0; i < tdim_; ++i)
    {
      // Set row to the global index
      row[i] = dofs[i];

      // Find position of the diagonal in the vectors
      // Not really sexy but we know that the diag index is in the vector,
      // so let us save a bound checking at each loop...
      // We cannot use a binary search as the indices might not be sorted.
      uint row_idx = dofs[i];
      uint diag_idx = 0;
      while (a_col_indices[i][diag_idx] != row_idx)
      {
        ++diag_idx;
      }

      // Scale the diagonal entry
      real diag_val = std::fabs(a[i][diag_idx]);

      // Zero the row
      uint nb_cols = a_col_indices[i].size();
      a_slip_row.resize(nb_cols, 0.0);
      As.set(&a_slip_row[0], 1, &row[i], nb_cols, &a_col_indices[i][0]);

      // Update the LHS diagonal
      As.set(&diag_val, 1, &row[i], 1, &row[i]);

      // Reset rhs for slip
      l_slip[i] = 0.0;
    }
  }

  //--- Apply local equation RHS to the copy of the matrix ---
  b.set(&l_slip[0], tdim_, &row[0]);

}

}
//-----------------------------------------------------------------------------

