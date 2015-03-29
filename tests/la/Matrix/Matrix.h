#ifdef HAVE_CHECK

#include <dolfin/la/Matrix.h>

#define MAT_SIZE 100

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_init_mat )
{
  Matrix A;
  A.init(MAT_SIZE, MAT_SIZE);

  fail_unless(A.size(0) == MAT_SIZE);
  fail_unless(A.size(1) == MAT_SIZE);

}END_TEST
//-----------------------------------------------------------------------------

Suite *test_suite_mat()
{

  TCase *tc;
  Suite *s;
  
  s = suite_create("Matrix");

  tc = tcase_create ("init");
  tcase_add_test (tc, test_init_mat);
  suite_add_tcase (s, tc);

  return s;
}


#endif
