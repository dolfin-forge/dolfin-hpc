#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include "ElementLibrary/ElementLibrary.h"

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
Suite *suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("elements");

  tc = tcase_create("ElementLibrary");
  tcase_add_test(tc, test_ElementLibrary);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  return s;
}
//-----------------------------------------------------------------------------
int main(void)
{

  int number_failed;
  Suite* s = suite();
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
  fprintf(stderr, "*** Check is required for dolfin/elements tests ***\n");
  return 0;
}

#endif
