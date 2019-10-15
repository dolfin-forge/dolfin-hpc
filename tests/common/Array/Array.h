#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Array.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Array )
  {
    {
      // default constructor, empty array of zero length
      Array<uint> A;
      ck_assert(A.ptr() == NULL);
      ck_assert(A.offset() == 0);
      ck_assert(A.stride() == 1);
    }
    //---
    {
      Array<uint> const B;
      ck_assert(B.ptr() == NULL);
      ck_assert(B.offset() == 0);
      ck_assert(B.stride() == 1);
    }
    //---
    {
      // array of given size
      Array<uint> C(2);
      C[0] = 1;
      C[1] = 2;
      ck_assert(C.ptr() != NULL);
      ck_assert(C.ptr() == &C[0]);
      ck_assert(C.size() == 2);
      ck_assert(*(C.ptr()) == C[0]);
      ck_assert(*(C.ptr() + 1) == C[1]);
    }
    //---
    {
      // array of given size with default value
      Array<uint> C(2,11);
      ck_assert(C.ptr() != NULL);
      ck_assert(C.ptr() == &C[0]);
      ck_assert(C.size() == 2);
      ck_assert(*(C.ptr()) == C[0]);
      ck_assert(*(C.ptr() + 1) == C[1]);
      ck_assert(C[0] == 11);
      ck_assert(C[1] == 11);
    }
    //---
    {
      // copy constructor + swap
      Array<uint> D(2);
      D[0] = 3;
      D[1] = 4;

      Array<uint> E(D);
      ck_assert(D.ptr() != E.ptr());
      ck_assert(D.size() == E.size());
      ck_assert(D.offset() == E.offset());
      ck_assert(D.stride() == E.stride());
      ck_assert(D[0] == E[0]);
      ck_assert(D[1] == E[1]);

      E[0] = 23;
      E[1] = 42;
      ck_assert(D.ptr() != E.ptr());
      ck_assert(D.size() == E.size());
      ck_assert(D.offset() == E.offset());
      ck_assert(D.stride() == E.stride());
      ck_assert(D[0] != E[0]);
      ck_assert(D[1] != E[1]);
      ck_assert(D[0] == 3);
      ck_assert(D[1] == 4);
      ck_assert(E[0] == 23);
      ck_assert(E[1] == 42);

      D.swap(E);
      ck_assert(D.ptr() != E.ptr());
      ck_assert(D.size() == E.size());
      ck_assert(D.offset() == E.offset());
      ck_assert(D.stride() == E.stride());
      ck_assert(D[0] != E[0]);
      ck_assert(D[1] != E[1]);
      ck_assert(D[0] == 23);
      ck_assert(D[1] == 42);
      ck_assert(E[0] == 3);
      ck_assert(E[1] == 4);
    }
    //---
    {
      //constructor from given range + operator= (copy and assign all elements)
      std::vector<uint> input(3,0);
      input[0] = 1;
      input[1] = 4;
      input[2] = 7;
      Array<uint> D(input.begin(), input.end());
      ck_assert(D[0] == 1);
      ck_assert(D[1] == 4);
      ck_assert(D[2] == 7);

      Array<uint> E(3,0);
      ck_assert(E[0] == 0);
      ck_assert(E[1] == 0);
      ck_assert(E[2] == 0);

      E = D;
      ck_assert(E[0] == 1);
      ck_assert(E[1] == 4);
      ck_assert(E[2] == 7);
      ck_assert(D[0] == 1);
      ck_assert(D[1] == 4);
      ck_assert(D[2] == 7);

      D = 2;
      ck_assert(D[0] == 2);
      ck_assert(D[1] == 2);
      ck_assert(D[2] == 2);
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
