#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/Simplex.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Simplex )
  {
    {
      Simplex<1> S;
      ck_assert(S.num_cells()    == 1);
      ck_assert(S.num_vertices() == 2);
      S.disp();
    }
    //---
    {
      Simplex<2> S;
      ck_assert(S.num_cells()    == 1);
      ck_assert(S.num_vertices() == 3);
      S.disp();
    }
    //---
    {
      Simplex<3> S;
      ck_assert(S.num_cells()    == 1);
      ck_assert(S.num_vertices() == 4);
      S.disp();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
