#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/mesh/Point.h>

#include <sstream>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<uint D>
void test_point()
{
  point<D> P;
  for(uint i = 0; i < D; ++i) { P[i] = 1 << D; }
  P.disp();
  P *= 2.0;
  P += 6.0;
  P -= 4.0;
  P /= 2.0;
  P.disp();
  P = 1.0;
  P.disp();
}
//-----------------------------------------------------------------------------
START_TEST( test_Point )
  {
    int init_failed = 0;
    begin("test_Point");
    //---
    {
      test_point<1>();
    }
    //---
    {
      test_point<2>();
    }
    //---
    {
      test_point<3>();
    }
    //---
    end();
    skip();
    ck_assert( init_failed == 0 );
  }END_TEST
//-----------------------------------------------------------------------------

#endif
