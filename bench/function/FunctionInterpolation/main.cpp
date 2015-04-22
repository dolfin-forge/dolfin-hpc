#include <dolfin.h>

#include "IthCoordinate.h"
#include "Coordinates.h"
#include "MPIRank.h"

#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/FunctionInterpolation.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);

  // Assume interpolation from M0 to M1
  // sqR is M0 and sqL is M1

  // Mesh generation
  //UnitSquare usqR(100,100);
  //File fsqR("square100R.xml");
  //fsqR << usqR;
  //UnitSquare usqL(100,100);
  //File fsqL("square100L.xml");
  //fsqL << usqL;
  //return 0;

  Mesh sqR("square100R.xml.gz");
  Mesh sqL("square100L.xml.gz");

  message("Scalar elements");
  {
    ufl::FiniteElement dg0(ufl::Family::DG, sqR.type(), 0);
    ufl::FiniteElement cg1(ufl::Family::CG, sqR.type(), 1);
    ufl::FiniteElement cg2(ufl::Family::CG, sqR.type(), 2);

    FiniteElementSpace dg0R(sqR, dg0);
    FiniteElementSpace cg1R(sqR, cg1);
    FiniteElementSpace cg2R(sqR, cg2);

    FiniteElementSpace dg0L(sqL, dg0);
    FiniteElementSpace cg2L(sqL, cg2);
    FiniteElementSpace cg1L(sqL, cg1);

    Mesh& M0 = sqR;
    Mesh& M1 = sqL;

    IthCoordinate X(M0, 0);

    FiniteElementSpace const& Vh0 = cg1R;
    Function F0(Vh0);
    Function F2(Vh0);

    FiniteElementSpace const& Vh1 = cg2L;
    Function F1(Vh1);
    Function F3(Vh1);

    message("Interpolate X to F0");
    FunctionInterpolation IX0(X, F0);
    IX0.compute();
    File f0("F0_S.pvd");
    f0 << F0;

    message("Interpolate F0 to F1");
    FunctionInterpolation I01(F0, F1);
    I01.compute();
    File f1("F1_S.pvd");
    f1 << F1;

    message("Interpolate F1 to F2");
    FunctionInterpolation I12(F1, F2);
    I12.compute();
    File f2("F2_S.pvd");
    f2 << F2;

    message("Interpolate X to F3");
    FunctionInterpolation IX3(X, F3);
    IX3.compute();
    File f3("F3_S.pvd");
    f3 << F3;

    real l20 = F0.vector().norm(dolfin::l2);
    F0.vector() -= F2.vector();
    F0.vector().apply();
    real l2err02 = F0.vector().norm(dolfin::l2);
    message("F0 - F2 l2 error = %16f", l2err02);

    real l21 = F1.vector().norm(dolfin::l2);
    F1.vector() -= F3.vector();
    F1.vector().apply();
    real l2err13 = F1.vector().norm(dolfin::l2);
    message("F1 - F3 l2 error = %16f", l2err13);
  }

  message("Vector elements");
  {
    uint const gdim = sqR.geometry().dim();
    ufl::VectorElement dg0d(ufl::Family::DG, sqR.type(), 0, gdim);
    ufl::VectorElement cg1d(ufl::Family::CG, sqR.type(), 1, gdim);
    ufl::VectorElement cg2d(ufl::Family::CG, sqR.type(), 2, gdim);

    FiniteElementSpace dg0dR(sqR, dg0d);
    FiniteElementSpace cg1dR(sqR, cg1d);
    FiniteElementSpace cg2dR(sqR, cg2d);

    FiniteElementSpace dg0dL(sqL, dg0d);
    FiniteElementSpace cg1dL(sqL, cg1d);
    FiniteElementSpace cg2dL(sqL, cg2d);

    Mesh& M0 = sqR;
    Mesh& M1 = sqL;

    Coordinates X(M0);

    FiniteElementSpace const& Vh0 = cg1dR;
    Function F0(Vh0);
    Function F2(Vh0);

    FiniteElementSpace const& Vh1 = cg2dL;
    Function F1(Vh1);
    Function F3(Vh1);

    message("Interpolate X to F0");
    FunctionInterpolation IX0(X, F0);
    IX0.compute();
    File f0("F0_V.pvd");
    f0 << F0;

    message("Interpolate F0 to F1");
    FunctionInterpolation I01(F0, F1);
    I01.compute();
    File f1("F1_V.pvd");
    f1 << F1;

    message("Interpolate F1 to F2");
    FunctionInterpolation I12(F1, F2);
    I12.compute();
    File f2("F2_V.pvd");
    f2 << F2;

    message("Interpolate X to F3");
    FunctionInterpolation IX3(X, F3);
    IX3.compute();
    File f3("F3_S.pvd");
    f3 << F3;

    real l20 = F0.vector().norm(dolfin::l2);
    F0.vector() -= F2.vector();
    F0.vector().apply();
    real l2err02 = F0.vector().norm(dolfin::l2);
    message("F0 - F2 l2 error = %16f", l2err02);

    real l21 = F1.vector().norm(dolfin::l2);
    F1.vector() -= F3.vector();
    F1.vector().apply();
    real l2err13 = F1.vector().norm(dolfin::l2);
    message("F1 - F3 l2 error = %16f", l2err13);
  }

  message("That's all folks !");

  return ret;
}
