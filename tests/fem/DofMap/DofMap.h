#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/DofMap.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_DofMap )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_DofMap");
  {
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
