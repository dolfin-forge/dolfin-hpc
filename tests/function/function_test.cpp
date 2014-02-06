#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>

#include <dolfin/function/Function.h>
#include <dolfin/function/ConstantFunction.h>
#include <dolfin/function/DiscreteFunction.h>
#include <dolfin/function/ExpressionFunction.h>
#include <dolfin/function/Expression.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/mesh/UnitSquare.h>

#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>

using dolfin::real;
using dolfin::message;
using dolfin::ConstantFunction;
using dolfin::DiscreteFunction;
using dolfin::ExpressionFunction;
using dolfin::Expression;
using dolfin::RealReference;

using ufl::Cell;
using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

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
START_TEST(test_init_discrete)
{
  int init_failed = 0;

  uint const deg_max = 2;
  std::vector<Family::Type> v;
  v.push_back(Family::DG);
  v.push_back(Family::CG);

  uint const N = 16;
  dolfin::UnitSquare mesh2d(N, N);

  for (std::vector<Family::Type>::const_iterator it = v.begin(); it != v.end();
      ++it)
  {
    Family f(*it);
    uint d_min = f.degree_min();
    uint d_max = std::min(f.degree_max(), std::max(d_min, deg_max));
    ufl::Domain::Set domains = f.domains();

    for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
        dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      Cell cell(dom);
      // Test just for UnitSquare
      if (cell.topological_dimension() == 2)
      {
        for (uint d = d_min; d <= d_max; ++d)
        {
          ufl::FiniteElement uflfem(*it, cell, d);
          DiscreteFunction F(mesh2d, uflfem.repr());
          F.disp();
        }
      }
    }
  }

  fail_unless(init_failed == 0);
}
END_TEST

//-----------------------------------------------------------------------------
START_TEST( test_init_expression )
{
  int init_failed = 0;

  uint const N = 16;
  dolfin::UnitSquare mesh2d(N, N);

  real rval = 1.0;
  RealReference expr(rval);
  ExpressionFunction F(mesh2d, expr);
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

  tcase_add_test(tc, test_init_constant);
  tcase_add_test(tc, test_init_discrete);
  tcase_add_test(tc, test_init_expression);

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
