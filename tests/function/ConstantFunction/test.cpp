#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>

#include <dolfin/function/ConstantFunction.h>
#include <dolfin/mesh/UnitSquare.h>

using dolfin::real;
using dolfin::message;
using dolfin::ConstantFunction;

#include <iostream>
#include <iomanip>

#ifdef HAVE_CHECK

#include <check.h>

int argc;
char * argv;

void setup()
{
}

void teardown()
{
}

//-----------------------------------------------------------------------------
START_TEST( test_init_constant )
{
  int init_failed = 0;

  uint const N = 16;
  dolfin::UnitSquare mesh2d(N, N);

  ConstantFunction F(mesh2d, 1.0);
  F.disp();

  fail_unless( init_failed == 0 );
}END_TEST

//-----------------------------------------------------------------------------

Suite *ufl_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("FUNCTION");
  tc = tcase_create("function");

  tcase_set_timeout(tc, 16);
  tcase_add_test(tc, test_init_constant);

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = ufl_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/function tests ***\n");
  return 0;
}

#endif
