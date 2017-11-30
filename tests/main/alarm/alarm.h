#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/main/alarm.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_alarm )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_alarm");
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
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
