#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Constant.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Constant )
  {
    Constant r;
    Constant const r0(0.0);
    Constant const r1(1.0);
    Constant const r2(2.0);

    // Default
    ck_assert(r == r0);

    r = 1.0;
    ck_assert(r == r1);

    r += 1.0;
    ck_assert(r == r2);

    r -= 1.0;
    ck_assert(r == r1);

    r *= 2.0;
    ck_assert(r == r2);

    r /= 2.0;
    ck_assert(r == r1);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
