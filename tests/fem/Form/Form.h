#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>
#include <dolfin/log/log.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/Functional.h>
#include <dolfin/fem/LinearForm.h>
#include <dolfin/fem/BilinearForm.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Functional )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Functional");
  {
  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_LinearForm )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_LinearForm");
  {
  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_BilinearForm )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_BilinearForm");
  {
  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
