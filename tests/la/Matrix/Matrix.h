#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/la/Matrix.h>

#define MAT_SIZE 100

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_init_mat )
{
  Matrix A;
  A.init(MAT_SIZE, MAT_SIZE);

  ck_assert(A.size(0) == MAT_SIZE);
  ck_assert(A.size(1) == MAT_SIZE);

}END_TEST
//-----------------------------------------------------------------------------

DOLFIN_SUITE_BEGIN(test_suite_mat, "Matrix")
{
  DOLFIN_TCASE_CREATE ("init");
  DOLFIN_TCASE_ADD(test_init_mat);
}
DOLFIN_SUITE_END

#endif
