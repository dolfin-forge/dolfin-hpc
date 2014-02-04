#include <dolfin/config/dolfin_config.h>

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>

using dolfin::FiniteElement;
using ufl::Cell;
using ufl::Domain;
using ufl::Family;
using ufl::Object;
using ufl::Space;

#include <iostream>

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

START_TEST( test_init_function_from_signature )
  {
    int init_failed = 0;

    uint const deg_max = 2;
    std::vector<Family::Type> v;
    v.push_back(Family::DG);
    v.push_back(Family::CG);

    for (std::vector<Family::Type>::const_iterator it = v.begin();
        it != v.end(); ++it)
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
        for (uint d = d_min; d <= d_max; ++d)
        {
          ufl::FiniteElement fem(*it, cell, d);
          std::cout << fem.repr() << std::endl;
        }
      }
    }

    fail_unless( init_failed == 0 );
  }END_TEST

Suite *fem_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("FEM");
  tc = tcase_create("fem");

  tcase_add_test(tc, test_init_function_from_signature);

  suite_add_tcase(s, tc);
  tcase_add_checked_fixture(tc, setup, teardown);

  return s;
}

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

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/fem tests ***\n");
  return 0;
}

#endif
