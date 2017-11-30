#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/main/SubSystemsManager.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_SubSystemsManager )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_SubSystemsManager");
  {
    SubSystemsManager::timer().set_limit(10);
    while (true)
    {
      if (SubSystemsManager::timer().state())
      {
        message("alarm::limit is true");
        break;
      }
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
