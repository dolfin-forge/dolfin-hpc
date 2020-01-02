#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/Functional.h>
#include <dolfin/fem/LinearForm.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/mesh/UnitSquare.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class T>
void test_nil()
{
  UnitSquare mesh( 5, 5 );
  CoefficientMap coefs;
  Nil<T>(mesh, coefs);
}

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Functional )
  {
    test_nil<Functional>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_LinearForm )
  {
    test_nil<LinearForm>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_BilinearForm )
  {
    test_nil<BilinearForm>();
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
