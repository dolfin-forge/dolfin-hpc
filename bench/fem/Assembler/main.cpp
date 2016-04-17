#include <dolfin.h>

#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/UFCCellIterator.h>

// 2D
#include "PoissonP1.h"
#include "PoissonP2.h"
#include "PoissonP3.h"
#include "StabStokes2D.h"
#include "THStokes2D.h"

// 3D
#include "Elasticity3D.h"
#include "NSEMomentum3D.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void test_numbering(Form& a)
{
  message("test_numbering");
  ufc::dofmap * ufc_dofmap = a.create_dofmap(0);
  DofNumbering * dn = DofNumbering::create(a.mesh(), *ufc_dofmap);
  dn->build();
  uint * dofs = new uint[ufc_dofmap->local_dimension()];
  for (UFCCellIterator cell(a.mesh()); !cell.end(); ++cell)
  {
    dn->tabulate_dofs(dofs, cell);
  }
  delete[] dofs;
  delete dn;
  delete ufc_dofmap;
}
//-----------------------------------------------------------------------------
void test_sparsity(Form& a)
{
  message("test_sparsity");
  UFC ufc(a);
  SparsityPattern sp;
  sp.init(ufc.form.rank(), ufc.global_dimensions, ufc.local_sizes);
  for (CellIterator cell(a.mesh()); !cell.end(); ++cell)
  {
    ufc.cell.update(*cell);
    for (uint i = 0; i < ufc.form.rank(); ++i)
    {
      a.dofmaps()[i].tabulate_dofs(ufc.dofs[i], ufc.cell);
    }
    sp.insert(ufc.local_dimensions, ufc.dofs);
  }
  sp.apply();
}
//-----------------------------------------------------------------------------
void assemble_a(Form& a, uint n)
{
  Assembler assembler;
  Matrix matA;
  tic();
  for (uint i = 0; i < n; ++i)
  {
    assembler.assemble(matA, a, true);
  }
  uint const N = a.dofmaps()[0].global_dimension();
  message("Assembly done with %u degrees of freedom", N);
  tocd();
}

int main(int argc, char** argv)
{
  dolfin_init(argc, argv);
  uint const n = 10;
  logm.verbose(1);
  logm.file();

  {
    Mesh mesh("../../data/meshes/squareN100R.xml.gz");

    //-------------------------------------------------------------------------
    {
      message("PoissonP1");
      PoissonP1BilinearForm a(mesh);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }

    //-------------------------------------------------------------------------
    {
      message("PoissonP2");
      PoissonP2BilinearForm a(mesh);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }

    //-------------------------------------------------------------------------
    {
      message("PoissonP3");
      PoissonP3BilinearForm a(mesh);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }

    //-------------------------------------------------------------------------
    {
      message("StabStokes2D");
      MeshSize h(mesh);
      StabStokes2DBilinearForm a(h);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }

    //-------------------------------------------------------------------------
    {
      message("THStokes2D");
      THStokes2DBilinearForm a(mesh);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }
  }

  {
    Mesh mesh("../../data/meshes/cubeN10R.xml.gz");
    //-------------------------------------------------------------------------
    {
      message("Elasticity3D");
      Elasticity3DBilinearForm a(mesh);
      test_numbering(a);
      test_sparsity(a);
      assemble_a(a, n);
    }

    //-------------------------------------------------------------------------
    {
      message("NSEMomentum3D");
      begin("");
      message("- Create functions.");
      Function um(mesh);
      Function delta1(mesh);
      Function delta2(mesh);
      Function k(mesh, 0.1);
      Function nu(mesh, 0.01);
      message("- Create form.");
      NSEMomentum3DBilinearForm a(um, delta1, delta2, k, nu);
      test_numbering(a);
      test_sparsity(a);
      message("- Create spaces.");
      FiniteElementSpace Sum(mesh, a, 2);
      um.init(Sum);
      FiniteElementSpace Sdelta1(mesh, a, 3);
      delta1.init(Sdelta1);
      FiniteElementSpace Sdelta2(mesh, a, 4);
      delta2.init(Sdelta2);
      message("- Assemble.");
      assemble_a(a, n);
      end();

    }
  }

  dolfin_finalize();
  return 0;
}

