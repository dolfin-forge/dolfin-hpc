#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/main/SubSystemsManager.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_SubSystemsManager )
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
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
