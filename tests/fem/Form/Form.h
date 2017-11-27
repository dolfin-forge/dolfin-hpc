#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/Test.h>
#include <dolfin/log/log.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/Functional.h>
#include <dolfin/fem/LinearForm.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/CoefficientMap.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class T>
void test_nil()
{
  Mesh mesh; CoefficientMap coefs; Nil<T>(mesh, coefs);
}

//-----------------------------------------------------------------------------
START_TEST( test_Functional )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Functional");
  {
    test_nil<Functional>();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_LinearForm )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_LinearForm");
  {
    test_nil<LinearForm>();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_BilinearForm )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_BilinearForm");
  {
    test_nil<BilinearForm>();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
