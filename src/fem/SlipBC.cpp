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
    mesh(mesh),
    sub_domains(NULL),
    sub_domain(0),
    sub_domains_local(false),
    user_sub_domain(&sub_domain),
    normal(mesh, VertexNormal::none),
    boundary(NULL),
    cell_map(NULL),
    vertex_map(NULL),
    As(NULL),
    row_block(NULL),
    zero_block(NULL),
    a1_indices_array(NULL),
    tdim_(mesh.topology().dim())
{
  // Initialize sub domain markers
  init(sub_domain);

  sub_system = SubSystem(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(VertexNormal& normal, const SubDomain& sub_domain) :
    BoundaryCondition("Slip"),
    mesh(normal.mesh()),
    sub_domains(NULL),
    sub_domain(0),
    sub_domains_local(false),
    user_sub_domain(&sub_domain),
    normal(normal),
    boundary(NULL),
    cell_map(NULL),
    vertex_map(NULL),
    As(NULL),
    row_block(NULL),
    zero_block(NULL),
    a1_indices_array(NULL),
    tdim_(mesh.topology().dim())
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
    normal(mesh, VertexNormal::none),
    boundary(NULL),
    cell_map(NULL),
    vertex_map(NULL),
    As(NULL),
    N_local(0),
    N_offset(0),
    row_block(NULL),
    zero_block(NULL),
    a1_indices_array(NULL),
    tdim_(mesh.topology().dim())
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, const SubDomain& sub_domain, const SubSystem& sub_system) :
    BoundaryCondition("Slip"),
    mesh(mesh),
    sub_domains(NULL),
    sub_domain(0),
    sub_domains_local(false),
    sub_system(sub_system),
    user_sub_domain(&sub_domain),
    normal(mesh, VertexNormal::none),
    boundary(NULL),
    cell_map(NULL),
    vertex_map(NULL),
    As(NULL),
    N_local(0),
    N_offset(0),
    row_block(NULL),
    zero_block(NULL),
    a1_indices_array(NULL),
    tdim_(mesh.topology().dim())
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
    normal(normal),
    boundary(NULL),
    cell_map(NULL),
    vertex_map(NULL),
    As(NULL),
    N_local(0),
    N_offset(0),
    row_block(NULL),
    zero_block(NULL),
    a1_indices_array(NULL),
    tdim_(mesh.topology().dim())
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
VertexNormal& SlipBC::normals()
{
  return normal;
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
  dolfin::error("Not implemented:\n",
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

  UFC ufc(form.form(), mesh, form.dofmaps());
  if (boundary == 0)
  {
    boundary = new BoundaryMesh(mesh);
    if (boundary->numCells())
    {
      cell_map = boundary->data().meshFunction("cell map");
      vertex_map = boundary->data().meshFunction("vertex map");
    }
  }
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
        ufc.update(*c, mesh.distdata());
        (form.dofmaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, c->index());

        for (uint j = 0; j < (form.dofmaps())[0].local_dimension(); j++)
          off_proc_rows.insert(ufc.dofs[0][j]);
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
          if (cvi[ci] == node)
            break;

        ufc.update(cell, mesh.distdata());
        (form.dofmaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, cell.index());

        // Get components of the vector-valued function at the current node.
        for (uint i = 0; i < gdim; i++, ci += cdim)
          nodes.push_back(ufc.dofs[0][ci]);

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
void SlipBC::init(const SubDomain& sub_domain)
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
    Array<uint>& dofs)
{

  // Node type defines the number of discriminated surface at the node.
  // Therefore it is the number of constrained directions up to the topological
  // dimension
  uint const n_type = std::max(normal.vertex_type.get(node), tdim_);

  for (uint i = 0; i < tdim_; ++i)
  {
    // Copy rows from the stiffnes matrix into lhs
    A.getrow(dofs[i], a_col_indices[i], a[i]);

    // Copy rhs to local
    l[i] = b[dofs[i]];

    // Fill node basis and determine maximum component
    basis_[i][0] = normal.basis()[i]->get(dofs[0]);
    maxcomp[i] = 0;
    for (uint j = 1; j < tdim_; ++j)
    {
      if ((basis_[i][j] = normal.basis()[i]->get(dofs[j]))
          > basis_[i][maxcomp[i]])
      {
        maxcomp[i] = j;
      }
    }
  }

//  uint * rows =

//
//  // Constrain directions
//  if (n_type >= tdim_)
//  {
//    for (uint i = 0; i < tdim_; ++i)
//    {
//      As.set(zero_block, 1, &dofs[i], static_cast<uint>(a_ncols[0]),
//          a1_indices_array);
//      Aset(As, dofs[i], dofs[i], 1.0);
//      bset(b, dofs[i], 0.0);
//    }
//
//    // Set rows of matrix to the identity matrix
//    for (uint i = 0; i < a_ncols[0]; i++)
//      row_block[i] = a[0][i] * tau2[0] + a[1][i] * tau2[1] + a[2][i] * tau2[2];
//
//    As.set(row_block, 1, &r3, static_cast<uint>(a_ncols[0]), a1_indices_array);
//    bset(b, r3, l[0] * tau2[0] + l[1] * tau2[1] + l[2] * tau2[2]);
//
//    for (uint i = 0; i < a_ncols[0]; i++)
//      row_block[i] = a[0][i] * tau1[0] + a[1][i] * tau1[1] + a[2][i] * tau1[2];
//
//    As.set(row_block, 1, &row_d1, static_cast<uint>(a_ncols[0]),
//        a1_indices_array);
//    bset(b, row_d1, l[0] * tau1[0] + l[1] * tau1[1] + l[2] * tau1[2]);
//
//    // Set second row of slip BC (tau..) to the row r3
//    for (uint i = 0; i < a_ncols[0]; i++)
//      row_block[i] = a[0][i] * tau2[0] + a[1][i] * tau2[1] + a[2][i] * tau2[2];
//
//    As.set(row_block, 1, &row_d2, static_cast<uint>(a_ncols[0]),
//        a1_indices_array);
//    bset(b, row_d2, l[0] * tau2[0] + l[1] * tau2[1] + l[2] * tau2[2]);
//
//  }

  // Apply local equation to the copy of the stiffness
  for (uint i = 0; i < tdim_; ++i)
  {
    As.setrow(dofs[i], a_col_indices[i], a[i]);
  }
}

}
//-----------------------------------------------------------------------------

