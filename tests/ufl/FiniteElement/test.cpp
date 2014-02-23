#include <dolfin/config/dolfin_config.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>

using namespace dolfin;

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
START_TEST( test_init )
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
        ufl::Cell cell(dom);
        for (uint d = d_min; d <= d_max; ++d)
        {
          begin("Creating UFLFiniteElement:");
          ufl::FiniteElement uflfem(*it, cell, d);
          message(uflfem.repr());
          uflfem.display();
          end();
          skip();

          begin("Creating UFLFiniteElement from factory function:");
          ufl::FiniteElementBase * factuflfem = ufl::FiniteElementBase::create(uflfem.repr());
          message(factuflfem->repr());
          factuflfem->display();
          if(factuflfem->repr() != uflfem.repr())
          {
              init_failed += 1;
          }
          delete factuflfem;
          end();
          skip();

        }
      }
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

Suite *ufl_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("UFL");
  tc = tcase_create("ufl");

  tcase_add_test(tc, test_init );

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
  fprintf(stderr, "*** Check is required for dolfin/ufl tests ***\n");
  return 0;
}

#endif
