
#ifndef DOLFIN_CHECK_H_
#define DOLFIN_CHECK_H_

#include <dolfin/main/init.h>

#include <stdio.h>

#ifdef HAVE_CHECK
#include <check.h>
#else
using Suite = void;
#endif

#include <dolfin/common/Test.h>
#include <dolfin/common/maybe_unused.h>

namespace dolfin
{

extern "C" typedef void(*CheckVoidFunctionPtr)(void);
#if (CHECK_MINOR_VERSION > 12)
extern "C" typedef const TTest * CheckIntFunctionPtr;
#else
extern "C" typedef void(*CheckIntFunctionPtr)(int);
#endif

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
  Suite *s = suite_create(_name); TCase *tc = nullptr;

#define DOLFIN_SUITE_END \
  return s; \
}

#define DOLFIN_SUITE_RUN(_name, _suite) \
  dolfin::Check::run_suite(_name, _suite);

#define DOLFIN_CHECK_SUITE(_name, _suite_function) \
int main() \
{ \
  return DOLFIN_SUITE_RUN(_name, _suite_function()); \
}

#if (CHECK_MINOR_VERSION > 12)
#define DOLFIN_START_TEST(_name) \
START_TEST( _name ) \
{\
  int init_failed = 0; \
  Test T;
#else
#define DOLFIN_START_TEST(_name) \
START_TEST( _name ) \
  int init_failed = 0; \
  Test T;
#endif

#if (CHECK_MINOR_VERSION > 12)
#define DOLFIN_END_TEST \
  ck_assert( init_failed == 0 ); \
} \
END_TEST
#else
#define DOLFIN_END_TEST \
  ck_assert( init_failed == 0 ); \
END_TEST
#endif

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

    MAYBE_UNUSED(name);
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
