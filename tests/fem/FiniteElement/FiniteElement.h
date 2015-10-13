#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/FiniteElement.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_FiniteElement )
{
  int init_failed = 0;
  Test::begin("test_FiniteElement");
  //---

  //---
  Test::end();
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
