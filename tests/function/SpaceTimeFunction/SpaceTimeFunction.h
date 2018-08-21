#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/SpaceTimeFunction.h>
#include <dolfin/main/PE.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_SpaceTimeFunction )
  {
    pushd(testresu("function", "SpaceTimeFunction"));
    Simplex<1> M;
    pushd("0");
    {
      SpaceTimeFunction S("U", FiniteElementSpace(M, DG(M.type(), 0)));
      real const dt = 0.1;
      for (Time t(0., 1.); t.is_valid(); t += dt)
      {
        S(t) = real(t);
        ck_assert(S.min() == S.max());
        ck_assert(S.min() == real(t));
        S.save();
        //message("(%e, %e)", S.min(), S.max());
      }
      for (Time t(0., 1.); t.is_valid(); t += dt)
      {
        S(t);
        S.eval();
        if (PE::rank() == 0)
        {
          message("%e = (%e, %e)", real(t), S.min(), S.max());
        }
        ck_assert(S.min() == S.max());
        ck_assert(S.min() == real(t));
      }
    }
    popd();
    pushd("1");
    {
      real const dt = 0.1;
      {
        SpaceTimeFunction S("U", FiniteElementSpace(M, DG(M.type(), 0)));
        for (Time t(0., 1.); t.is_valid(); t += dt)
        {
          S(t) = real(t);
          ck_assert(S.min() == S.max());
          ck_assert(S.min() == real(t));
          S.save();
          //message("(%e, %e)", S.min(), S.max());
        }
      }
      {
        SpaceTimeFunction S("U", FiniteElementSpace(M, DG(M.type(), 0)));
        for (Time t(0., 1.); t.is_valid(); t += dt)
        {
          S(t);
          S.eval();
          if (PE::rank() == 0)
          {
            message("%e = (%e, %e)", real(t), S.min(), S.max());
          }
          ck_assert(S.min() == S.max());
          ck_assert(S.min() == real(t));
        }
      }
    }
    popd();
    popd();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
