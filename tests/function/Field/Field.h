#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Field.h>
#include <dolfin/mesh/Simplex.h>
#include <dolfin/ufl/UFLFiniteElement.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Field )
{
  int init_failed = 0;
  Test T;
  //---
  uint const D = 2;       // Field depth
  uint const N = 10;      // Number of steps
  //---
  T.begin("test_Field");
  {
    Simplex<1> m;
    FiniteElementSpace Vh(m, scalar<DG>(m.type(), 0));
    Field U(Vh);
    U.allocate(D);
    for (uint d = 0; d <= D; ++d)
    {
      U[d] = real(d);
    }
    U.disp();
    //
    message("Field shift");
    tic();
    for (uint i = 0; i < N; ++i)
    {
      U.shift();
      U.disp();
    }
    tocd();
    //
    message("Field vector copy");
    tic();
    for (uint n = 0; n < N; ++n)
    {
      for (uint d = D; d > 0; --d)
      {
        U[d] = U[d - 1];
        U[d].sync();
      }
    }
    tocd();
    //
    message("Field 1-shift");
    tic();
    for (uint d = 0; d <= D; ++d)
    {
      U[d] = real(d);
    }
    U.disp();
    for (uint i = 0; i < 8; ++i)
    {
      U << 1;
      U.disp();
    }
    tocd();
    //
    message("Field 2-shift");
    tic();
    for (uint d = 0; d <= D; ++d)
    {
      U[d] = real(d);
    }
    U.disp();
    for (uint i = 0; i < 4; ++i)
    {
      U << 2;
      U.disp();
    }
    tocd();
  }
  T.end();
  //---
  T.begin("Field with timestamp");
  {
    Simplex<1> m;
    FiniteElementSpace Vh(m, scalar<DG>(m.type(), 0));
    Field U(Vh);
    U.allocate(D);
    for (uint d = 0; d <= D; ++d)
    {
      U[d] = real(d);
    }
    U.disp();
    //
    message("Field shift");
    tic();
    Time t(0.0, 1.0);
    U[0](t);
    for (uint i = 1; i <= N; ++i)
    {
      U.shift();
      t.clock() = i*t.measure()/real(N);
      U[0](t);
      U.disp();
    }
    tocd();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
