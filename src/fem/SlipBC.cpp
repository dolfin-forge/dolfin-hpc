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
SlipBC::SlipBC(Mesh& mesh, SubDomain& sub_domain) :
        BoundaryCondition("Slip"),
        mesh(mesh),
        sub_domains(NULL),
        sub_domain(0),
        sub_domains_local(false),
        user_sub_domain(&sub_domain),
        node_normal(mesh),
        As(NULL),
        row_block(NULL),
        zero_block(NULL),
        a1_indices_array(NULL),
        boundary(NULL),
        cell_map(NULL),
        vertex_map(NULL)
{
  // Fill permutation matrix
  permutation_matrix_[0][0] = 0;
  permutation_matrix_[0][1] = 1;
  permutation_matrix_[0][2] = 2;
  permutation_matrix_[1][0] = 1;
  permutation_matrix_[1][1] = 0;
  permutation_matrix_[1][2] = 2;
  permutation_matrix_[2][0] = 2;
  permutation_matrix_[2][1] = 1;
  permutation_matrix_[2][2] = 0;

  // Initialize sub domain markers
  init(sub_domain);

  sub_system = SubSystem(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain& sub_domain, NodeNormal& node_normal) :
        BoundaryCondition("Slip"),
        mesh(mesh),
        sub_domains(NULL),
        sub_domain(0),
        sub_domains_local(false),
        user_sub_domain(&sub_domain),
        node_normal(node_normal),
        As(NULL),
        row_block(NULL),
        zero_block(NULL),
        a1_indices_array(NULL),
        boundary(NULL),
        cell_map(NULL),
        vertex_map(NULL)
{
  // Initialize sub domain markers
  init(sub_domain);

  sub_system = SubSystem(0);
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain) :
        BoundaryCondition("Slip"),
        mesh(sub_domains.mesh()),
        sub_domains(&sub_domains),
        sub_domain(sub_domain),
        sub_domains_local(false),
        sub_system(),
        user_sub_domain(NULL),
        node_normal(mesh),
        nzm(0),
        As(NULL),
        N_local(0),
        N_offset(0),
        row_block(NULL),
        zero_block(NULL),
        a1_indices_array(NULL),
        boundary(NULL),
        cell_map(NULL),
        vertex_map(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain& sub_domain, const SubSystem& sub_system) :
        BoundaryCondition("Slip"),
        mesh(mesh),
        sub_domains(NULL),
        sub_domain(0),
        sub_domains_local(false),
        sub_system(sub_system),
        user_sub_domain(&sub_domain),
        node_normal(mesh),
        nzm(0),
        As(NULL),
        N_local(0),
        N_offset(0),
        row_block(NULL),
        zero_block(NULL),
        a1_indices_array(NULL),
        boundary(NULL),
        cell_map(NULL),
        vertex_map(NULL)
{
  // Set sub domain markers
  init(sub_domain);
  //Function g(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
    const SubSystem& sub_system) :
        BoundaryCondition("Slip"),
        mesh(sub_domains.mesh()),
        sub_domains(&sub_domains),
        sub_domain(sub_domain),
        sub_domains_local(false),
        sub_system(sub_system),
        user_sub_domain(NULL),
        node_normal(node_normal),
        nzm(0),
        As(NULL),
        N_local(0),
        N_offset(0),
        row_block(NULL),
        zero_block(NULL),
        a1_indices_array(NULL),
        boundary(NULL),
        cell_map(NULL),
        vertex_map(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::~SlipBC()
{
  // Delete sub domain markers if created locally
  if (sub_domains_local)
    delete sub_domains;

  if (As)
    delete As;

  if (a1_indices_array)
    delete[] a1_indices_array;
  if (row_block)
    delete[] row_block;
  if (zero_block)
    delete[] zero_block;
  if (boundary)
    delete boundary;
}
//-----------------------------------------------------------------------------
void
SlipBC::apply(GenericMatrix& A, GenericVector& b, const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void
SlipBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map,
    const ufc::form& ufc_form)
{
  dolfin::error("Not implemented:\n",
      "void apply(GenericMatrix& A, GenericVector& b,\n",
      "\tDofMap const& dof_map, const ufc::form& form)");
}
//-----------------------------------------------------------------------------
void
SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
    const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void
SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
    const DofMap& dof_map, const ufc::form& ufc_form)
{
  dolfin::error("Not implemented:\n",
      "void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,\n",
      "\tDofMap const& dof_map, const ufc::form& form)");
}

//-----------------------------------------------------------------------------
void
SlipBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map,
    const Form& form)
{

  if (MPI::processNumber() == 0)
  {
    dolfin_set("output destination", "terminal");
  }
  message("Applying SlipBC boundary conditions to linear system.");
  dolfin_set("output destination", "silent");

  UFC ufc(form.form(), mesh, form.dofMaps());

  /// Create boundary mesh
  if (boundary == 0)
  {
    boundary = new BoundaryMesh(mesh);
    if (boundary->numCells())
    {
      cell_map = boundary->data().meshFunction("cell map");
      vertex_map = boundary->data().meshFunction("vertex map");
    }
  }

  /// Create matrix to hold a copy of the current linear system
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
      /// Initialize the list of DoFs not on process
      std::map<uint, uint> mapping;
      for (CellIterator c(mesh); !c.end(); ++c)
      {
        ufc.update(*c, mesh.distdata());
        (form.dofMaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, c->index());

        for (uint j = 0; j < (form.dofMaps())[0].local_dimension(); j++)
        {
          off_proc_rows.insert(ufc.dofs[0][j]);
        }
      }

      b.init_ghosted(off_proc_rows.size(), off_proc_rows, mapping);
    }

    row_block = new real[A.size(0)];
    zero_block = new real[A.size(0)];
    a1_indices_array = new uint[A.size(0)];
  }
  // Copy global stiffness matrix into temporary one
  *(As->instance()) = A;

  Array<uint> nodes;
  uint gdim = mesh.geometry().dim();
  uint cdim = mesh.type().numVertices(mesh.topology().dim());

  uint count = 0;
  if (boundary->numCells())
  {
    /// Loop on boundary cells
    for (VertexIterator v(*boundary); !v.end(); ++v)
    {
      Vertex vertex(mesh, vertex_map->get(*v));

      // Skip facets not inside the sub domain
      if ((*sub_domains)(vertex) != sub_domain)
      {
        continue;
      }

      uint node = vertex.index();
      if (!mesh.distdata().is_ghost(node, 0) || MPI::numProcesses() == 1)
      {
        Cell cell(mesh, (vertex.entities(gdim))[0]);

        uint *cvi = cell.entities(0);
        uint ci = 0;
        for (ci = 0; ci < cell.numEntities(0); ci++)
        {
          if (cvi[ci] == node)
          {
            break;
          }
        }

        ufc.update(cell, mesh.distdata());
        (form.dofMaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, cell.index());

        // Get components of the vector-valued function at the current node.
        //FIXME: Assumes that the value dimension is the geometrical dimension
        for (uint i = 0; i < gdim; i++, ci += cdim)
        {
          nodes.push_back(ufc.dofs[0][ci]);
        }
        applySlipBC((Matrix&) A, *As, (Vector&) b, mesh, node, nodes);
        count++;
        nodes.clear();
      }
    }
  }

  // Apply changes in the temporary matrix
  As->apply();

  // Apply changes in the stiffness matrix and load vector
  A = *(As->instance());
  b.apply();

}
//-----------------------------------------------------------------------------
void
SlipBC::init(SubDomain& sub_domain)
{
  // Create mesh function for sub domain markers on facets
  mesh.init(0);
  sub_domains = new MeshFunction<uint>(mesh, 0);
  sub_domains_local = true;

  // Mark everything as sub domain 1
  (*sub_domains) = 1;

  // Mark the sub domain as sub domain 0
  sub_domain.mark(*sub_domains, 0);
}

//-----------------------------------------------------------------------------
void
SlipBC::applySlipBC(Matrix& A, Matrix& As, Vector& b, Mesh& mesh, uint node,
    Array<uint>& nodes)
{

  // Now, we use node_type vector, which is defined for all nodes at the boundary:
  // node_type = 1: node at the surface; Apply slip BC as usual
  // node_type = 2: node at the edge; Apply slip BC in modified way using the second normal which is tau_1
  // node_type = 3: node at the corner; Apply no slip BC, meaning that u = 0 at this node

  // Note a[1] is just index(a[0]) + N

  // define type of node, in 3D: 1 - surface; 2 - edge; 3 - corner;
  //                      in 2D: 1 - surface; 2 - corner;

  uint const nsdim = mesh.topology().dim();
  uint const nb_nodes = nodes.size();

  Array<real> a[3];
  Array<uint> a_idx[3];
  uint ncols[3];
  std::fill_n(ncols, 3, 0);
  real l[3];

  // Get non-zero entries on row and get rhs
  for (uint d = 0; d < nsdim; ++d)
  {
    A.getrow(nodes[d], a_idx[d], a[d]);
    ncols[d] = a_idx[d].size();
    l[d] = b[nodes[d]];
  }

  // Attributes
  nzm = ncols[0];
  std::fill_n(row_block, nzm, 0.0);
  std::fill_n(zero_block, nzm, 0.0);

  // Copy non-zero entries column indices to a1_indices_array
  std::copy(a_idx[0].begin(), a_idx[0].end(), a1_indices_array);
  uint n_type = node_normal.node_type.get(node);

  // Apply no-slip for the node on the corner in 3D and 2D
  if (n_type >= nsdim)
  {
    for (uint d = 0; d < nsdim; ++d)
    {
      As.set(zero_block, 1, &nodes[d], static_cast<uint>(ncols[0]),
          a1_indices_array);
      Aset(As, nodes[d], nodes[d], 1.0);
      bset(b, nodes[d], 0.0);
    }
  }
  else
  {
    std::fill_n(n_, 3, 0.0);
    std::fill_n(tau_1_, 3, 0.0);
    std::fill_n(tau_2_, 3, 0.0);

    real maxn = 0.0;
    uint maxrow = 0;
    for (uint d = 0; d < nsdim; ++d)
    {
      n_[d] = node_normal.normal[d].get(node);
      tau_1_[d] = node_normal.tau_1[d].get(node);
      tau_2_[d] = node_normal.tau_2[d].get(node);

      // Find maximum absolute value of the normal components
      if (fabs(n_[d]) > maxn)
      {
        maxn = fabs(n_[d]);
        maxrow = d;
      }
    }

    // row[0], row[1], row[2] are rows which are corresponds to the boundary point
    // find the maximum component of the normal and put it to the diagonal
    uint const (&row)[3] = permutation_matrix_[maxrow];

    As.set(zero_block, 1, &row[0], static_cast<uint>(ncols[0]),
        a1_indices_array);

    Aset(As, row[0], nodes[0], n_[0]);
    Aset(As, row[0], nodes[1], n_[1]);
    if (nsdim == 3)
      Aset(As, row[0], nodes[2], n_[2]);

    bset(b, row[0], 0.0);

    // Set [nx,ny,nz] on row with maxn and b[maxn] = 0

    // set the row using tangent vector in 2D
    /*
     if (nsdim == 2) {
     rows[0] = row[1];

     // Set rows of matrix to the identity matrix
     for(uint i = 0; i < ncols[0]; i++)
     row_block[i] = (a[0][i] * t1 + a[1][i] * t2);

     As.set(row_block, 1, rows, static_cast<uint>(ncols[0]), a1_indices_array);

     bset(b, row[1], (l[0] * t1 + l[1] * t2));
     }
     */

    // if node is on the edge in 3D case:
    if (n_type == 2 && nsdim == 3)
    {

      // The variables below are used in rearranging rows
      real nn[3];
      real tt[3];

      // Count number of zero entries in the normal vector
      int non0 = 0;
      for (uint d = 0; d < nsdim; ++d)
      {
        if (fabs(n_[d]) < DOLFIN_EPS)
          non0++;
      }

      if (non0 > 1)
      {
        // nn = n_ and tt = tau_1_
        nn[0] = n_[0];
        nn[1] = n_[1];
        nn[2] = n_[2];
        tt[0] = tau_1_[0];
        tt[1] = tau_1_[1];
        tt[2] = tau_1_[2];
      }
      else
      {
        // tt = n_ and nn = n_
        tt[0] = n_[0];
        tt[1] = n_[1];
        tt[2] = n_[2];
        nn[0] = tau_1_[0];
        nn[1] = tau_1_[1];
        nn[2] = tau_1_[2];
      }

      // find maximum of the normal components and put them to the diagonal
      real maxn = nn[0];
      uint maxrow = 0;
      for (uint d = 1; d < nsdim; ++d)
      {
        // Find maximum absolute value of the normal components
        if (fabs(nn[d]) > maxn)
        {
          maxn = fabs(n_[d]);
          maxrow = d;
        }
      }

      // row[0], row[1], row[2] are rows which are corresponds to the boundary point
      // find the maximum component of the normal and put it to the diagonal
      uint row[3];
      row[0]= permutation_matrix_[maxrow][0];
      row[1]= permutation_matrix_[maxrow][1];
      row[2]= permutation_matrix_[maxrow][2];

      As.set(zero_block, 1, &row[0], static_cast<uint>(ncols[0]),
          a1_indices_array);

      Aset(As, row[0], nodes[0], nn[0]);
      Aset(As, row[0], nodes[1], nn[1]);
      Aset(As, row[0], nodes[2], nn[2]);

      maxn = 0.0;
      maxn = fmax(fabs(tt[0]), fabs(tt[1]));
      maxn = fmax(maxn, fabs(tt[2]));

      // Start rearanging rows
      if (fabs(fabs(tt[0]) - maxn) < DOLFIN_EPS)
      { // tt[0] is a largest component
        row[1] = nodes[0];
        row[2] = (row[0] == nodes[1] ? nodes[2] : nodes[1]);
        if (row[1] == row[0])
        {
          real mm = fmax(fabs(tt[1]), fabs(tt[2]));
          if (fabs(fabs(tt[1]) - mm) < DOLFIN_EPS)
          {
            row[1] = nodes[1];
            row[2] = nodes[2];
          }
          else
          {
            row[1] = nodes[2];
            row[2] = nodes[1];
          }
        }
      }
      else if (fabs(fabs(tt[1]) - maxn) < DOLFIN_EPS)
      { // tt[1] is a largest component
        row[1] = nodes[1];
        row[2] = (row[0] == nodes[0] ? nodes[2] : nodes[0]);
        if (row[1] == row[0])
        {
          real mm = fmax(fabs(tt[0]), fabs(tt[2]));
          if (fabs(fabs(tt[0]) - mm) < DOLFIN_EPS)
          {
            row[1] = nodes[0];
            row[2] = nodes[2];
          }
          else
          {
            row[1] = nodes[2];
            row[2] = nodes[0];
          }
        }
      }
      else
      { // tt[2] is a largest component
        row[1] = nodes[2];
        row[2] = (row[0] == nodes[0] ? nodes[1] : nodes[0]);
        if (row[1] == row[0])
        {
          real mm = fmax(fabs(tt[0]), fabs(tt[1]));
          if (fabs(fabs(tt[0]) - mm) < DOLFIN_EPS)
          {
            row[1] = nodes[0];
            row[2] = nodes[1];
          }
          else
          {
            row[1] = nodes[1];
            row[2] = nodes[0];
          }
        }
      }

      As.set(zero_block, 1, &row[1], static_cast<uint>(ncols[0]),
          a1_indices_array);

      Aset(As, row[1], nodes[0], tt[0]);
      Aset(As, row[1], nodes[1], tt[1]);
      Aset(As, row[1], nodes[2], tt[2]);
      bset(b, row[1], 0.0);

      // Set second row of slip BC (tau..) to the row row[2]

      // Set rows of matrix to the identity matrix
      for (uint i = 0; i < ncols[0]; i++)
        row_block[i] = a[0][i] * tau_2_[0] + a[1][i] * tau_2_[1]
            + a[2][i] * tau_2_[2];

      //As.set(row_block, rows, 1, a1_indices_array, ncols[0]);
      As.set(row_block, 1, &row[2], static_cast<uint>(ncols[0]),
          a1_indices_array);
      bset(b, row[2], l[0] * tau_2_[0] + l[1] * tau_2_[1] + l[2] * tau_2_[2]);

    }

    // The case where node is on the surface
    if (n_type == 1 && nsdim == 3)
    {
      int ind1 = 0;
      int ind2 = 0;

      for (uint i = 0; i < ncols[0]; i++)
      {
        if (row[1] == a_idx[0][i])
          ind1 = i;
        if (row[2] == a_idx[0][i])
          ind2 = i;
      }

      // find maximum of the elements of 2 vectors which are belong to diagonal
      real row_d11 = a[0][ind1] * tau_1_[0] + a[1][ind1] * tau_1_[1]
          + a[2][ind1] * tau_1_[2];
      real row_d21 = a[0][ind1] * tau_2_[0] + a[1][ind1] * tau_2_[1]
          + a[2][ind1] * tau_2_[2];
      real row_d12 = a[0][ind2] * tau_1_[0] + a[1][ind2] * tau_1_[1]
          + a[2][ind2] * tau_1_[2];
      real row_d22 = a[0][ind2] * tau_2_[0] + a[1][ind2] * tau_2_[1]
          + a[2][ind2] * tau_2_[2];

      real maxr = fmax(fabs(row_d11), fabs(row_d21));
      maxr = fmax(maxr, fabs(row_d12));
      maxr = fmax(maxr, fabs(row_d22));

      // define new rows according the above maximum:
      uint row_d1 = 0;
      uint row_d2 = 0;
      if (fabs(fabs(row_d11) - maxr) < DOLFIN_EPS
          || fabs(fabs(row_d22) - maxr) < DOLFIN_EPS)
      {
        row_d1 = row[1];
        row_d2 = row[2];
      }
      else
      {
        row_d1 = row[2];
        row_d2 = row[1];
      }

      for (uint i = 0; i < ncols[0]; i++)
      {
        row_block[i] = a[0][i] * tau_1_[0] + a[1][i] * tau_1_[1]
            + a[2][i] * tau_1_[2];
      }

      //As.set(row_block, rows, 1, a1_indices_array, ncols[0]);
      As.set(row_block, 1, &row_d1, static_cast<uint>(ncols[0]),
          a1_indices_array);
      bset(b, row_d1, l[0] * tau_1_[0] + l[1] * tau_1_[1] + l[2] * tau_1_[2]);

      // Set second row of slip BC (tau..) to the row row[2]
      for (uint i = 0; i < ncols[0]; i++)
      {
        row_block[i] = a[0][i] * tau_2_[0] + a[1][i] * tau_2_[1]
            + a[2][i] * tau_2_[2];
      }

      // Set row for u_tau2
      // A(tau2,j) = a[0][j] * tau_2_[0] + a[1][j] * tau_2_[1] + a[2][j] * tau_2_[2]
      // b(tau2)   = l[0] * tau_2_[0] + l[1] * tau_2_[1] + l[2] * tau_2_[2]

      //As.set(row_block, rows, 1, a1_indices_array, ncols[0]);
      As.set(row_block, 1, &row_d2, static_cast<uint>(ncols[0]),
          a1_indices_array);
      bset(b, row_d2, l[0] * tau_2_[0] + l[1] * tau_2_[1] + l[2] * tau_2_[2]);

    }
  }
  delete[] ncols;
}

}
//-----------------------------------------------------------------------------

