
#ifndef DOLFIN_CHECK_H_
#define DOLFIN_CHECK_H_

#include <dolfin/main/init.h>

#include <stdio.h>

#ifdef HAVE_CHECK
#include <check.h>
#else
typedef void Suite;
#endif

namespace dolfin
{

extern "C" typedef void(*CheckVoidFunctionPtr)(void);
extern "C" typedef void(*CheckIntFunctionPtr)(int);


//-----------------------------------------------------------------------------

#define DOLFIN_TCASE_CREATE(_name) \
  tc = tcase_create(_name); \
  suite_add_tcase(s, tc); \
  tcase_add_checked_fixture(tc, \
    (dolfin::CheckVoidFunctionPtr) dolfin::Check::setup, \
    (dolfin::CheckVoidFunctionPtr) dolfin::Check::teardown); \
  tcase_set_timeout(tc, 60);

#define DOLFIN_TCASE_ADD(_test) \
  tcase_add_test(tc, (dolfin::CheckIntFunctionPtr) _test);

#define DOLFIN_TCASE_TIMEOUT(_value) \
  tcase_set_timeout(tc, _value);

#define DOLFIN_SUITE_BEGIN(_suite_function, _name) \
Suite *_suite_function() \
{ \
  Suite *s = suite_create(_name); TCase *tc = NULL;

#define DOLFIN_SUITE_END \
  return s; \
}

#define DOLFIN_SUITE_RUN(_name, _suite) \
  dolfin::Check::run_suite(_name, _suite);

#define DOLFIN_CHECK_SUITE(_name, _suite_function) \
int main(int argc, char **argv) \
{ \
  return DOLFIN_SUITE_RUN(_name, _suite_function()); \
}

//-----------------------------------------------------------------------------

struct Check
{

  //---------------------------------------------------------------------------
  static void setup()
  {
  }

  //---------------------------------------------------------------------------
  static void teardown()
  {
  }

  //---------------------------------------------------------------------------
  static int run_suite(const char * name, Suite * s)
  {
#ifdef HAVE_CHECK

    int number_failed;
    SRunner* sr = srunner_create(s);
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? 0 : 1;

#else

    fprintf(stderr, "*** Check is required for running %s tests ***\n", name);
    return 1;

#endif
  }

  //---------------------------------------------------------------------------

};

} /* namespace dolfin */

#endif /* DOLFIN_CHECK_H_ */
