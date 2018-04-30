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
  P *= 2.0;
  P += 6.0;
  P -= 4.0;
  P /= 2.0;
  P = 1.0;
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Point )
  {
    test_point<1>();
    test_point<2>();
    test_point<3>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
