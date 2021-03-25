#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/la/Matrix.h>

#define MAT_SIZE 100

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_init_mat )
{
  Matrix A;
  A.init(MAT_SIZE, MAT_SIZE);

  ck_assert(A.size(0) == MAT_SIZE);
  ck_assert(A.size(1) == MAT_SIZE);

}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
#ifdef HAVE_TRILINOS
#include <dolfin/la/trilinos/TrilinosMatrix.h>
DOLFIN_START_TEST( test_trilinos_mat )
{
  trilinos::Matrix x1;
  trilinos::Matrix x2( MAT_SIZE, MAT_SIZE );
}DOLFIN_END_TEST
#else
DOLFIN_START_TEST( test_trilinos_mat )
{
}DOLFIN_END_TEST
#endif
//-----------------------------------------------------------------------------

#endif
