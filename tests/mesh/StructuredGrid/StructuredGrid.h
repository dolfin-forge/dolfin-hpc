#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/BoundingBox.h>
#include <dolfin/mesh/StructuredGrid.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_BoundingBox )
  {
    for (uint i = 0; i <= Space::MAX_DIMENSION; ++i)
    {
      BoundingBox bb(i);
      Point u(1.0, 1.0, 1.0);
      bb += u;
      bb -= u;
      real h = 2.0;
      bb *= h;
      Point d(0.5, 0.25, 0.125);
      bb *= d;
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_StructuredGrid_interval )
  {
    uint N = 8192;
    StructuredGrid g(IntervalCell(), N);
    dolfin_assert(g.num_cells() == N);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_StructuredGrid_triangle )
  {
    uint N = 128;
    StructuredGrid g(TriangleCell(), N);
    dolfin_assert(g.num_cells() == N*N*2);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_StructuredGrid_tetrahedron )
  {
    uint N = 32;
    StructuredGrid g(TetrahedronCell(), N);
    dolfin_assert(g.num_cells() == N*N*N*6);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_StructuredGrid_quadrilateral )
  {
    uint N = 128;
    StructuredGrid g(QuadrilateralCell(), N);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_StructuredGrid_hexahedron )
  {
    uint N = 32;
    StructuredGrid g(HexahedronCell(), N);
    dolfin_assert(g.num_cells() == N*N*N);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
