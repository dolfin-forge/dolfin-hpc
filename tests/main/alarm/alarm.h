#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/main/alarm.h>
#include <dolfin/log/log.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_alarm )
  {
    alarm s;
    bool b = s.set_limit(10);
    if (!b)
    {
      error("alarm : not set correctly");
    }
    while (true)
    {
      if (s.state())
      {
        message("alarm::limit is true");
        break;
      }
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
