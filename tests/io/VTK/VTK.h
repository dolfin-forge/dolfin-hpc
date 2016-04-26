#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/io/VTKFile.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_VTKMesh )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_VTKMesh");
  {
  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
