#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <FiniteElement/FiniteElement.h>
#include <FiniteElementSpace/FiniteElementSpace.h>
#include <BoundaryNormal/BoundaryNormal.h>

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
Suite* fem_suite()
{
  TCase* tc;
  Suite* s;

  s = suite_create("fem");

  tc = tcase_create("FiniteElement");
  tcase_add_test(tc, test_FiniteElement);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc, 60);

  tc = tcase_create("FiniteElementSpace");
  tcase_add_test(tc, test_FiniteElementSpace);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc, 60);

  tc = tcase_create("BoundaryNormal");
  tcase_add_test(tc, test_NodeNormal);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc, 60);

  return s;
}
//-----------------------------------------------------------------------------
int main(void)
{
  int number_failed;
  Suite* s = fem_suite();
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
  fprintf(stderr, "*** Check is required for dolfin/fem tests ***\n");
  return 0;
}

#endif
