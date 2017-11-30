#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DiscreteSpaces.h>
#include <dolfin/mesh/CellTypes.h>

using namespace dolfin;

template<class Element, class CellType, uint degree>
void test()
{
  Element E(CellType(), degree);
  E.disp();
}

//-----------------------------------------------------------------------------
START_TEST( test_DiscreteSpaces )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_DiscreteSpaces::CG");
  {
    test<CG, IntervalCell   , 1>();
    test<CG, TriangleCell   , 1>();
    test<CG, TetrahedronCell, 1>();
    test<CG, IntervalCell   , 2>();
    test<CG, TriangleCell   , 2>();
    test<CG, TetrahedronCell, 2>();
  }
  T.end();
  //---
  T.begin("test_DiscreteSpaces::DG");
  {
    test<DG, IntervalCell   , 1>();
    test<DG, TriangleCell   , 1>();
    test<DG, TetrahedronCell, 1>();
    test<DG, IntervalCell   , 2>();
    test<DG, TriangleCell   , 2>();
    test<DG, TetrahedronCell, 2>();
  }
  T.end();
  //---
  T.begin("test_DiscreteSpaces::vector<CG>");
  {
  }
  T.end();
  //---
  T.begin("test_DiscreteSpaces::vector<DG>");
  {
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
