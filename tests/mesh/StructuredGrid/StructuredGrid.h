#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/BoundingBox.h>
#include <dolfin/mesh/StructuredGrid.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/io/VTKFile.h>

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
  StructuredGrid g(IntervalCell(), N);
  dolfin_assert(g.numCells() == N);
  VTKFile vtk("StructuredGrid_interval.pvd");
  vtk << g;
  BinaryFile bin("StructuredGrid_interval.bin");
  bin << g;
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
  StructuredGrid g(TriangleCell(), N);
  dolfin_assert(g.numCells() == N*N*2);
  VTKFile vtk("StructuredGrid_triangle.pvd");
  vtk << g;
  BinaryFile bin("StructuredGrid_triangle.bin");
  bin << g;
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
  StructuredGrid g(TetrahedronCell(), N);
  dolfin_assert(g.numCells() == N*N*N*6);
  VTKFile vtk("StructuredGrid_tetrahedron.pvd");
  vtk << g;
  BinaryFile bin("StructuredGrid_tetrahedron.bin");
  bin << g;
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
  StructuredGrid g(QuadrilateralCell(), N);
  VTKFile vtk("StructuredGrid_quadrilateral.pvd");
  vtk << g;
  BinaryFile bin("StructuredGrid_quadrilateral.bin");
  bin << g;
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
  StructuredGrid g(HexahedronCell(), N);
  dolfin_assert(g.numCells() == N*N*N);
  VTKFile vtk("StructuredGrid_hexahedron.pvd");
  vtk << g;
  BinaryFile bin("StructuredGrid_hexahedron.bin");
  bin << g;
  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
