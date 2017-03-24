#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/function/Constant.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Constant )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Constant");
  {
    Constant r;
    Constant const r0(0.0);
    Constant const r1(1.0);
    Constant const r2(2.0);

    // Default
    ck_assert( r == r0 );

    r = 1.0;
    ck_assert( r == r1 );

    r += 1.0;
    ck_assert( r == r2 );

    r -= 1.0;
    ck_assert( r == r1 );

    r *= 2.0;
    ck_assert( r == r2 );

    r /= 2.0;
    ck_assert( r == r1 );
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
