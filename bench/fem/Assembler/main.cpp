#include <dolfin.h>

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
  message("Number of degrees of freedom: %d", N);
  tocd();
}

int main(int argc, char** argv)
{
  UnitSquare mesh2D(128, 128);
  UnitCube mesh3D(32, 32, 32);
  uint n = 1;

  //---------------------------------------------------------------------------
  {
    message("PoissonP1");
    PoissonP1BilinearForm a(mesh2D);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("PoissonP2");
    PoissonP2BilinearForm a(mesh2D);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("PoissonP3");
    PoissonP3BilinearForm a(mesh2D);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("StabStokes2D");
    MeshSize h(mesh2D);
    StabStokes2DBilinearForm a(h);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("THStokes2D");
    THStokes2DBilinearForm a(mesh2D);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("Elasticity3D");
    Elasticity3DBilinearForm a(mesh3D);
    assemble_a(a, n);
  }

  //---------------------------------------------------------------------------
  {
    message("NSEMomentum3D");
    begin("");
    message("- Create functions.");
    Function um(mesh3D);
    Function delta1(mesh3D);
    Function delta2(mesh3D);
    Function k(mesh3D, 0.1);
    Function nu(mesh3D, 0.01);
    message("- Create form.");
    NSEMomentum3DBilinearForm a(um, delta1, delta2, k, nu);
    message("- Create spaces.");
    FiniteElementSpace Sum(mesh3D, a, 2);
    um.init(Sum);
    FiniteElementSpace Sdelta1(mesh3D, a, 3);
    delta1.init(Sdelta1);
    FiniteElementSpace Sdelta2(mesh3D, a, 4);
    delta2.init(Sdelta2);
    message("- Assemble.");
    assemble_a(a, n);
    end();

  }

  return 0;
}

