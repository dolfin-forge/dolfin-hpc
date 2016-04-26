#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include "Constant/Constant.h"
#include "Expression/Expression.h"
#include "Real/Real.h"
#include "Value/Value.h"
#include "UFCFunction/UFCFunction.h"

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

  s = suite_create("function");

  tc = tcase_create("Constant");
  tcase_add_test(tc, test_Constant);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("Expression");
  tcase_add_test(tc, test_Expression);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("Value");
  tcase_add_test(tc, test_Value);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("Real");
  tcase_add_test(tc, test_Real);
  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);
  tcase_set_timeout(tc,60);

  tc = tcase_create("UFCFunction");
  tcase_add_test(tc, test_UFCFunction);
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
  fprintf(stderr, "*** Check is required for dolfin/function tests ***\n");
  return 0;
}

#endif
