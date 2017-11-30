#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include "SubDomain/SubDomain.h"
#include "Point/Point.h"
#include "CellTypes/CellTypes.h"
#include "UnitMeshes/UnitMeshes.h"
#include "MeshFunction/MeshFunction.h"
#include "MeshData/MeshData.h"
#include "Simplex/Simplex.h"
#include "StructuredGrid/StructuredGrid.h"
#include "VertexNormal/VertexNormal.h"
#include "algorithm/algorithm.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "mesh")
{
  DOLFIN_TCASE_CREATE("SubDomain");
  DOLFIN_TCASE_ADD(test_SubDomain);

  DOLFIN_TCASE_CREATE("Point");
  DOLFIN_TCASE_ADD(test_Point);

  DOLFIN_TCASE_CREATE("CellTypes");
  DOLFIN_TCASE_ADD(test_PointCell);
  DOLFIN_TCASE_ADD(test_IntervalCell);
  DOLFIN_TCASE_ADD(test_TriangleCell);
  DOLFIN_TCASE_ADD(test_TetrahedronCell);
  DOLFIN_TCASE_ADD(test_QuadrilateralCell);
  DOLFIN_TCASE_ADD(test_HexahedronCell);

  DOLFIN_TCASE_CREATE("Simplex");
  DOLFIN_TCASE_ADD(test_Simplex);

  DOLFIN_TCASE_CREATE("UnitMeshes");
  DOLFIN_TCASE_ADD(test_UnitInterval);
  DOLFIN_TCASE_ADD(test_UnitSquare);
  DOLFIN_TCASE_ADD(test_UnitCube);
  DOLFIN_TCASE_ADD(test_Box);
  DOLFIN_TCASE_ADD(test_UnitDisk);

  DOLFIN_TCASE_CREATE("StructuredGrid");
  DOLFIN_TCASE_ADD(test_BoundingBox);
  DOLFIN_TCASE_ADD(test_StructuredGrid_interval);
  DOLFIN_TCASE_ADD(test_StructuredGrid_triangle);
  DOLFIN_TCASE_ADD(test_StructuredGrid_tetrahedron);
  DOLFIN_TCASE_ADD(test_StructuredGrid_quadrilateral);
  DOLFIN_TCASE_ADD(test_StructuredGrid_hexahedron);

  DOLFIN_TCASE_CREATE("MeshFunction");
  DOLFIN_TCASE_ADD(test_MeshFunction );

  DOLFIN_TCASE_CREATE("MeshData");
  DOLFIN_TCASE_ADD(test_MeshData_add );
  DOLFIN_TCASE_ADD(test_MeshData_insert );
  DOLFIN_TCASE_ADD(test_MeshData_iterator );

  DOLFIN_TCASE_CREATE("VertexNormal");
  DOLFIN_TCASE_ADD(test_VertexNormal);

  DOLFIN_TCASE_CREATE("algorithm");
  DOLFIN_TCASE_ADD(test_algorithm );
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/mesh", suite)
//-----------------------------------------------------------------------------

#endif
