#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>

#include <dolfin/function/Constant.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/UnitSquare.h>

using dolfin::real;
using dolfin::message;
using dolfin::Constant;
using dolfin::Point;

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

  real const ref0 = 1.0;
  Constant F(mesh2d, ref0);
  real value[1] = { 0.0 };
  Point x;
  F.eval( value, &x[0] );
  dolfin_assert(value[0] == ref0);

  real const ref1 = 2.0 * ref0;
  F = ref1;
  F.eval( value, &x[0] );
  dolfin_assert(value[0] == ref1);


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
