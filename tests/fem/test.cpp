#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/init.h>

#include <iostream>
#include <iomanip>

#ifdef HAVE_CHECK

#include <check.h>

#include <BoundaryNormal/BoundaryNormal.h>
#include <FiniteElement/FiniteElement.h>
#include <FiniteElementSpace/FiniteElementSpace.h>

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

Suite * test_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("FEM");

  tc = tcase_create("SubDomain");
  tcase_add_test(tc, test_nodenormal_create);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("FiniteElement");
  tcase_add_test(tc, test_init_element);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("FiniteElementSpace");
  tcase_add_test(tc, test_init_element_space);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = test_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/fem tests ***\n");
  return 0;
}

#endif
