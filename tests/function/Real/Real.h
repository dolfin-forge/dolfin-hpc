#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Real.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Real1 )
  {
    Real<> r;
    r.disp();
    r = 0.0;
    r.disp();
    ck_assert(r[0] == 0.0);
    r += 1.0;
    r.disp();
    ck_assert(r[0] == 1.0);
    r *= 2.0;
    r.disp();
    ck_assert(r[0] == 2.0);
    r /= 4.0;
    r.disp();
    ck_assert(r[0] == 0.5);
    r -= 0.25;
    r.disp();
    ck_assert(r[0] == 0.25);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Real2 )
  {
    Real<2> r;
    r.disp();
    r[0] = 0.0;
    r[1] = 1.0;
    r.disp();
    ck_assert(r[0] == 0.0);
    ck_assert(r[1] == 1.0);
    r += 1.0;
    r.disp();
    ck_assert(r[0] == 1.0);
    ck_assert(r[1] == 2.0);
    r *= 2.0;
    r.disp();
    ck_assert(r[0] == 2.0);
    ck_assert(r[1] == 4.0);
    r /= 4.0;
    r.disp();
    ck_assert(r[0] == 0.5);
    ck_assert(r[1] == 1.0);
    r -= 0.25;
    r.disp();
    ck_assert(r[0] == 0.25);
    ck_assert(r[1] == 0.75);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Realt )
  {
    Real<1> r(2.0);
    Time t(0.0, 1.0);
    real const k = t.measure() / 10;
    for (real s = t.begin(); t.is_valid(k/10); t.step(k))
    {
      (void) s;
      r(t);
      ck_assert(r[0] == 2.0);
      t.disp();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
