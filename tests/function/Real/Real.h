#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Real.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Real )
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
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
