#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/XMLFile.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_XMLMesh )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_XMLMesh");
  {
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
