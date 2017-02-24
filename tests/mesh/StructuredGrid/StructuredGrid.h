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
  begin("test_BoundingBox");
  //---
  for (dolfin::uint i = 0; i <= Space::MAX_DIMENSION; ++i)
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
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_interval )
{
  int init_failed = 0;
  begin("test_StructuredGrid_interval");
  //---
  dolfin::uint N = 8192;
  StructuredGrid g(IntervalCell(), N);
  dolfin_assert(g.num_cells() == N);
  //---
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_triangle )
{
  int init_failed = 0;
  begin("test_StructuredGrid_triangle");
  //---
  dolfin::uint N = 128;
  StructuredGrid g(TriangleCell(), N);
  dolfin_assert(g.num_cells() == N*N*2);
  //---
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_tetrahedron )
{
  int init_failed = 0;
  begin("test_StructuredGrid_tetrahedron");
  //---
  dolfin::uint N = 32;
  StructuredGrid g(TetrahedronCell(), N);
  dolfin_assert(g.num_cells() == N*N*N*6);
  //---
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_quadrilateral )
{
  int init_failed = 0;
  begin("test_StructuredGrid_quadrilateral");
  //---
  dolfin::uint N = 128;
  StructuredGrid g(QuadrilateralCell(), N);
  //---
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_StructuredGrid_hexahedron )
{
  int init_failed = 0;
  begin("test_StructuredGrid_hexahedron");
  //---
  dolfin::uint N = 32;
  StructuredGrid g(HexahedronCell(), N);
  dolfin_assert(g.num_cells() == N*N*N);
  //---
  end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
