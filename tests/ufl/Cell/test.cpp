#include <dolfin/config/dolfin_config.h>

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLDomain.h>

using ufl::Cell;
using ufl::Domain;

#include <iostream>
#include <iomanip>

using namespace dolfin;

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

    Domain::Set domains;
    domains.insert(domains.begin(), Domain::interval);
    domains.insert(domains.begin(), Domain::triangle);
    domains.insert(domains.begin(), Domain::tetrahedron);

    for (ufl::Domain::Set::const_iterator dom_it = domains.begin();
         dom_it != domains.end(); ++dom_it)
    {
      Domain dom(*dom_it);
      Cell cell(dom);
      cell.display();
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
