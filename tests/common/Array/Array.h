#include <dolfin/common/Check.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/common/Array.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Array )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Array::ptr");
  {
    Array<uint> A;
    ck_assert(A.ptr() == NULL);

    Array<uint> const B;
    ck_assert(B.ptr() == NULL);

    Array<uint> C(2); C[0] = 1; C[1] = 2;
    ck_assert(C.ptr() != NULL);
    ck_assert(C.ptr() == &C[0]);
    ck_assert(C.size() == 2);
    ck_assert(*(C.ptr()) == C[0]);
    ck_assert(*(C.ptr() + 1) == C[1]);
    uint const * c = C.ptr();
    for (uint i = 0; i < 2; ++i) message("%u", *(c++));

    Array<uint> D(2); D[0] = 3; D[1] = 4;
    ck_assert(D.ptr() != NULL);
    ck_assert(D.ptr() == &D[0]);
    ck_assert(D.size() == 2);
    ck_assert(*(D.ptr()) == D[0]);
    ck_assert(*(D.ptr() + 1) == D[1]);
    uint const * d = D.ptr();
    for (uint i = 0; i < 2; ++i) message("%u", *(d++));
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
