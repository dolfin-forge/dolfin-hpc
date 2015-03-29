#include <dolfin/config/dolfin_config.h>

#include <iostream>
#include <iomanip>

#ifdef HAVE_CHECK

#include <check.h>

#include "SubDomain/SubDomain.h"
#include "UnitMeshes/UnitMeshes.h"
#include "VertexNormal/VertexNormal.h"


int argc;
char **argv;

void setup()
{
  //  dolfin_init(argc, argv);
}

void teardown()
{
  //  dolfin_finalize();
}

Suite *mesh_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("Mesh");

  tc = tcase_create("SubDomain");
  tcase_add_test(tc, test_geometric_create);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("UnitMeshes");
  tcase_add_test(tc, test_create_interval);
  tcase_add_test(tc, test_create_square);
  tcase_add_test(tc, test_create_cube);
  tcase_add_test(tc, test_create_box);
  tcase_add_test(tc, test_create_disk);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("VertexNormal");
  tcase_add_test(tc, test_init_weight_none );
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  return s;
}

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

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/mesh tests ***\n");
  return 0;
}

#endif
