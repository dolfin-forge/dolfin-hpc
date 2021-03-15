#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Constant.h>
#include <dolfin/function/Analytic.h>
#include <dolfin/mesh/Simplex.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Analytic )
  {
    Simplex<1> M;
    Analytic<Constant> A(M);

    // Default
    //    ck_assert(real(evaluant(A)) == 0.0);

    // Modified value
    //    evaluant(A) = 1.0;
    //    ck_assert(real(evaluant(A)) == 1.0);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
