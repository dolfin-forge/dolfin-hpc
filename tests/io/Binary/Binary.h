#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/io/BinaryFile.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_BinaryMesh )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_BinaryMesh");
  {
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
