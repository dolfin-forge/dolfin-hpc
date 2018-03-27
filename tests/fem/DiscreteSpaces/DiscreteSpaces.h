#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DiscreteSpaces.h>
#include <dolfin/mesh/CellTypes.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class Element, class CellType, uint degree>
void test()
{
  Element E(CellType(), degree);
}

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_DiscreteSpaces )
  {
    {
      test<CG, IntervalCell   , 1>();
      test<CG, TriangleCell   , 1>();
      test<CG, TetrahedronCell, 1>();
      test<CG, IntervalCell   , 2>();
      test<CG, TriangleCell   , 2>();
      test<CG, TetrahedronCell, 2>();
    }
    //---
    {
      test<DG, IntervalCell   , 1>();
      test<DG, TriangleCell   , 1>();
      test<DG, TetrahedronCell, 1>();
      test<DG, IntervalCell   , 2>();
      test<DG, TriangleCell   , 2>();
      test<DG, TetrahedronCell, 2>();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
