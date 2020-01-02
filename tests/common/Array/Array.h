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
    //---
    {
      /// Assign to all elements in the array
      Array<uint> a(5, 42);
      a = 23;
      ck_assert(a[0] == 23);
      ck_assert(a[1] == 23);
      ck_assert(a[2] == 23);
      ck_assert(a[3] == 23);
      ck_assert(a[4] == 23);
      a = 42;
      ck_assert(a[0] == 42);
      ck_assert(a[1] == 42);
      ck_assert(a[2] == 42);
      ck_assert(a[3] == 42);
      ck_assert(a[4] == 42);
    }
    //---
    {
      /// Assignement operator
      Array<uint> a(5, 42);
      Array<uint> b = a;
      ck_assert(b.size() == 5);
      ck_assert(b.ptr() != a.ptr());
      ck_assert(b.offset() == a.offset());
      ck_assert(b.stride() == a.stride());
      ck_assert(a[0] == b[0]);
      ck_assert(a[1] == b[1]);
      ck_assert(a[2] == b[2]);
      ck_assert(a[3] == b[3]);
      ck_assert(a[4] == b[4]);
      Array<uint> c(7, 23);
      ck_assert(c.size() == 7);
      c = a;
      ck_assert(c.size() == 5);
      ck_assert(c.ptr() != a.ptr());
      ck_assert(c.offset() == a.offset());
      ck_assert(c.stride() == a.stride());
      ck_assert(a[0] == c[0]);
      ck_assert(a[1] == c[1]);
      ck_assert(a[2] == c[2]);
      ck_assert(a[3] == c[3]);
      ck_assert(a[4] == c[4]);
    }
    //---
    {
      /// append()
      Array<uint> a(2, 42);
      Array<uint> b(3, 23);
      ck_assert(a.size() == 2);
      a.append( b.begin(), b.end() );
      ck_assert(a.size() == 5);
      ck_assert(a.ptr() != b.ptr());
      ck_assert(a[0] == 42);
      ck_assert(a[1] == 42);
      ck_assert(a[2] == 23);
      ck_assert(a[3] == 23);
      ck_assert(a[4] == 23);

      ck_assert(b.size() == 3);
      bool threw_exception = false;
      try {
        b.append( a.end(), a.begin() );
      } catch ( ... )
      {
        threw_exception = true;
      }
      ck_assert(threw_exception == true);
      ck_assert(b.size() == 3);
    }
    //---
    // {
    //   /// Assignement operator

    // }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
