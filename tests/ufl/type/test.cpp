#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/types.h>

using dolfin::uint;
using dolfin::real;

#include <dolfin/ufl/UFLtype.h>

using ufl::type;

#include <iostream>
#include <iomanip>
#include <vector>

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
START_TEST( test_init_uint )
  {
    int init_failed = 0;

    for (uint d = 0; d < 4; ++d)
    {
      type<uint> t(d);
      t.display();
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_real )
  {
    int init_failed = 0;

    for (uint d = 0; d < 4; ++d)
    {
      real val = 1.1 * (real) d;
      type<real> t(val);
      t.display();
    }

    fail_unless( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_init_string )
  {
    int init_failed = 0;

    std::vector<std::string> vals;
    vals.push_back("'interval'");
    vals.push_back("'triangle'");
    vals.push_back("'tetrahedron'");

    for (std::vector<std::string>::const_iterator it = vals.begin();
        it != vals.end(); ++it)
    {
      type<std::string> t(*it);
      t.display();
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

  tcase_add_test(tc, test_init_uint);
  tcase_add_test(tc, test_init_real);
  tcase_add_test(tc, test_init_string);

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
