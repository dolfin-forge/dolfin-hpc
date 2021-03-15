#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Operators.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/celltypes/HexahedronCell.h>
#include <dolfin/mesh/celltypes/IntervalCell.h>
#include <dolfin/mesh/celltypes/QuadrilateralCell.h>
#include <dolfin/mesh/celltypes/TetrahedronCell.h>
#include <dolfin/mesh/celltypes/TriangleCell.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template < class CellType >
void check_cell_ops()
{
  CellType cellt;
  Mesh     cellm( cellt, EuclideanSpace( cellt.dim() ) );
  cellt.create_reference_cell( cellm );
  MeshValues< real, Cell > M( cellm );

  Circumradius< Cell > O0;
  M << O0;
  Diameter< Cell > O1;
  M << O1;
  Volume< Cell > O2;
  M << O2;
}

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Operators )
{
  check_cell_ops< IntervalCell >();
  check_cell_ops< TriangleCell >();
  check_cell_ops< TetrahedronCell >();
  check_cell_ops< QuadrilateralCell >();
  check_cell_ops< HexahedronCell >();
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
