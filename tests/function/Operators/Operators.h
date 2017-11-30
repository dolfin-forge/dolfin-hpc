#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Operators.h>
#include <dolfin/mesh/PointCell.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/QuadrilateralCell.h>
#include <dolfin/mesh/HexahedronCell.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class CellType>
void check_cell_ops()
{
  CellType cellt;
  Mesh     cellm;
  cellt.create_reference_cell(cellm);
  MeshValues<real, Cell> M(cellm);

  Circumradius<Cell> O0; M << O0;
  Diameter<Cell>     O1; M << O1;
  Volume<Cell>       O2; M << O2;
}

//-----------------------------------------------------------------------------
START_TEST( test_Operators )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Operators");
  {
    check_cell_ops<IntervalCell>();
    check_cell_ops<TriangleCell>();
    check_cell_ops<TetrahedronCell>();
    check_cell_ops<QuadrilateralCell>();
    check_cell_ops<HexahedronCell>();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
