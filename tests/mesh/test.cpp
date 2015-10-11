#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include "SubDomain/SubDomain.h"
#include "CellTypes/CellTypes.h"
#include "UnitMeshes/UnitMeshes.h"
#include "StructuredGrid/StructuredGrid.h"
#include "VertexNormal/VertexNormal.h"

#include <check.h>

//-----------------------------------------------------------------------------
void setup()
{
}
//-----------------------------------------------------------------------------
void teardown()
{
}
//-----------------------------------------------------------------------------
Suite *mesh_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("mesh");

  tc = tcase_create("SubDomain");
  tcase_add_test(tc, test_SubDomain);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("CellTypes");
  tcase_add_test(tc, test_PointCell);
  tcase_add_test(tc, test_IntervalCell);
  tcase_add_test(tc, test_TriangleCell);
  tcase_add_test(tc, test_TetrahedronCell);
  tcase_add_test(tc, test_QuadrilateralCell);
  tcase_add_test(tc, test_HexahedronCell);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("UnitMeshes");
  tcase_add_test(tc, test_UnitInterval);
  tcase_add_test(tc, test_UnitSquare);
  tcase_add_test(tc, test_UnitCube);
  tcase_add_test(tc, test_Box);
  tcase_add_test(tc, test_UnitDisk);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("StructuredGrid");
  tcase_add_test(tc, test_BoundingBox);
  tcase_add_test(tc, test_StructuredGrid_interval);
  tcase_add_test(tc, test_StructuredGrid_triangle);
  tcase_add_test(tc, test_StructuredGrid_tetrahedron);
  tcase_add_test(tc, test_StructuredGrid_quadrilateral);
  tcase_add_test(tc, test_StructuredGrid_hexahedron);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("VertexNormal");
  tcase_add_test(tc, test_VertexNormal );
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  return s;
}
//-----------------------------------------------------------------------------
int main(void)
{

  int number_failed;
  Suite* s = mesh_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;

}
//-----------------------------------------------------------------------------
#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/mesh tests ***\n");
  return 0;
}

#endif
