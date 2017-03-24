#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/function/Real.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Real )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Value : E01");
  {
    Real<> r;
    r.disp();
    r = 0.0;
    r.disp();
    r += 1.0;
    r.disp();
    r *= 2.0;
    r.disp();
    r /= 4.0;
    r.disp();
    r -= 0.25;
    r.disp();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
