#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/EntityKey.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<uint D>
void check_map()
{
  std::map<EntityKey, uint> M;
  EntityKey k0(D);
  EntityKey k1(D);
  EntityKey k2(D);
  for (uint i = 0; i < D; ++i)
  {
    k0.indices[i] = i;
    k1.indices[i] = i * 2;
    k2.indices[i] = i + 1;
  }
  k0.disp();
  k1.disp();
  k2.disp();
  ck_assert(k0 == k0);
  ck_assert(k1 == k1);
  ck_assert(k2 == k2);
  ck_assert(k0 <  k1);
  ck_assert(k0 <= k1);
  ck_assert(k1 <  k2);
  ck_assert(k1 <= k2);
  ck_assert(k2 >  k1);
  ck_assert(k2 >= k1);
  ck_assert(k1 >  k0);
  ck_assert(k1 >= k0);

  M[k0] = 0;
  ck_assert(M.count(k0));
  ck_assert(M.find(k0) != M.end());
  ck_assert(M.size() == 1);
  M[k1] = 1;
  ck_assert(M.count(k1));
  ck_assert(M.find(k1) != M.end());
  ck_assert(M.size() == 2);
  M[k2] = 2;
  ck_assert(M.count(k2));
  ck_assert(M.find(k2) != M.end());
  ck_assert(M.size() == 3);

  for (std::map<EntityKey, uint>::const_iterator it = M.begin(); it != M.end();
       ++it)
  {
    it->first.disp();
  }
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_EntityKey )
  {
    check_map<2>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
