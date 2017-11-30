#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Expression.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Expression )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Expression");
  {
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
