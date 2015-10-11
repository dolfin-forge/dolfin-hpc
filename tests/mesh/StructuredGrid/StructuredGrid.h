#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/BoundingBox.h>
#include <dolfin/mesh/StructuredGrid.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_BoundingBox )
{
  int init_failed = 0;
  Test::begin("test_BoundingBox");
  //---
  for (uint i = 0; i <= EuclideanSpace::MAX_DIMENSION; ++i)
  {
    BoundingBox bb(i);
    bb.disp();
    Point u(1.0, 1.0, 1.0);
    bb += u;
    bb.disp();
    bb -= u;
    real h = 2.0;
    bb *= h;
    bb.disp();
    Point d(0.5, 0.25, 0.125);
    bb *= d;
    bb.disp();
  }
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_interval )
{
  int init_failed = 0;
  Test::begin("test_StructuredGrid_interval");
  //---
  uint N = 8192;
  StructuredGrid g1(IntervalCell(), N);
  dolfin_assert(g1.numCells() == N);
  VTKFile vtk1("StructuredGrid_interval.pvd");
  vtk1 << g1;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_triangle )
{
  int init_failed = 0;
  Test::begin("test_StructuredGrid_triangle");
  //---
  uint N = 128;
  StructuredGrid g2(TriangleCell(), N);
  dolfin_assert(g2.numCells() == N*N*2);
  VTKFile vtk2("StructuredGrid_triangle.pvd");
  vtk2 << g2;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_tetrahedron )
{
  int init_failed = 0;
  Test::begin("test_StructuredGrid_tetrahedron");
  //---
  uint N = 32;
  StructuredGrid g3(TetrahedronCell(), N);
  dolfin_assert(g3.numCells() == N*N*N*6);
  VTKFile vtk3("StructuredGrid_tetrahedron.pvd");
  vtk3 << g3;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_quadrilateral )
{
  int init_failed = 0;
  Test::begin("test_StructuredGrid_quadrilateral");
  //---
  uint N = 128;
  StructuredGrid g4(QuadrilateralCell(), N);
  VTKFile vtk4("StructuredGrid_quadrilateral.pvd");
  vtk4 << g4;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_hexahedron )
{
  int init_failed = 0;
  Test::begin("test_StructuredGrid_hexahedron");
  //---
  uint N = 32;
  StructuredGrid g6(HexahedronCell(), N);
  dolfin_assert(g6.numCells() == N*N*N);
  VTKFile vtk6("StructuredGrid_hexahedron.pvd");
  vtk6 << g6;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
