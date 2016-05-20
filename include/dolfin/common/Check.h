
#ifndef DOLFIN_CHECK_H_
#define DOLFIN_CHECK_H_

#include <dolfin/common/types.h>
#include <dolfin/main/init.h>

#include <stdio.h>

#ifdef HAVE_CHECK
#include <check.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
static void check_setup()
{
}

//-----------------------------------------------------------------------------
static void check_teardown()
{
}

//-----------------------------------------------------------------------------
static int check_run_suite(Suite * s)
{
#ifdef HAVE_CHECK
  dolfin::dolfin_init(0, NULL);
  int number_failed;
  SRunner* sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  dolfin::dolfin_finalize();
  return (number_failed == 0) ? 0 : 1;

#else

  fprintf(stderr, "*** Check is required for running tests ***\n");
  return 1;

#endif
}

//-----------------------------------------------------------------------------
#ifdef HAVE_CHECK
START_TEST( test_dummy )
{
  int init_failed = 0;
  //---
  fprintf(stdout, "*** Check dummy test ***\n");
  //---
  fail_unless( init_failed == 0 );
}END_TEST
#endif
//-----------------------------------------------------------------------------

}

#endif /* DOLFIN_CHECK_H_ */
