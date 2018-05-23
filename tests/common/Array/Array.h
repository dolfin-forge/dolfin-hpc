#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Array.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Array )
  {
    {
      Array<uint> A;
      ck_assert(A.ptr() == NULL);
    }
    //---
    {
      Array<uint> const B;
      ck_assert(B.ptr() == NULL);
    }
    //---
    {
      Array<uint> C(2);
      C[0] = 1;
      C[1] = 2;
      ck_assert(C.ptr() != NULL);
      ck_assert(C.ptr() == &C[0]);
      ck_assert(C.size() == 2);
      ck_assert(*(C.ptr()) == C[0]);
      ck_assert(*(C.ptr() + 1) == C[1]);
      C.dump();
    }
    //---
    {
      Array<uint> D(2);
      D[0] = 3;
      D[1] = 4;
      ck_assert(D.ptr() != NULL);
      ck_assert(D.ptr() == &D[0]);
      ck_assert(D.size() == 2);
      ck_assert(*(D.ptr()) == D[0]);
      ck_assert(*(D.ptr() + 1) == D[1]);
      D.dump();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
