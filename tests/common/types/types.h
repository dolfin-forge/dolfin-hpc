#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/types.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_types )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_types::intersection");
  {
    {
      _set<uint> s0; s0.insert(0);
      _set<uint> s1; s1.insert(1);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.empty());
    }

    {
      _set<uint> s0; s0.insert(1);
      _set<uint> s1; s1.insert(1);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.size() == 1 && in.count(1));
    }

    {
      _set<uint> s0; s0.insert(1); s0.insert(2); s0.insert(3);
      _set<uint> s1; s1.insert(0); s1.insert(1); s1.insert(2);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.size() == 2 && in.count(1) && in.count(2));
    }

    {
      _set<uint> s0; s0.insert(0); s0.insert(1); s0.insert(2);
      _set<uint> s1; s1.insert(1); s1.insert(2); s1.insert(3);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.size() == 2 && in.count(1) && in.count(2));
    }

    {
      _set<uint> s0; s0.insert(1); s0.insert(2); s0.insert(3); s0.insert(4);
      _set<uint> s1; s1.insert(0); s1.insert(1); s1.insert(4);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.size() == 2 && in.count(1) && in.count(4));
    }

    {
      _set<uint> s0; s0.insert(0); s0.insert(1); s0.insert(4);
      _set<uint> s1; s1.insert(1); s1.insert(2); s1.insert(3); s1.insert(4);
      _set<uint> in;
      intersection(s0, s1, in);
      ck_assert(in.size() == 2 && in.count(1) && in.count(4));
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
