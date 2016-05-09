#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/UnitInterval.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Function )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Function");
  {
    UnitInterval mesh(42);





  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
