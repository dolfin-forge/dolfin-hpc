#ifdef HAVE_CHECK

#include <dolfin/la/Vector.h>

#define VEC_SIZE 100

using namespace dolfin;
//-----------------------------------------------------------------------------
START_TEST( test_init_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  fail_unless(x.size() == VEC_SIZE);
  fail_unless(x.local_size() == VEC_SIZE);

}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_zero_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  x = 1.0;
  x.zero();
  fail_unless(x.max() == 0);
  
}END_TEST
//-----------------------------------------------------------------------------

Suite *test_suite_vec()
{

  TCase *tc;
  Suite *s;
  
  s = suite_create("Vector");

  tc = tcase_create ("init");
  tcase_add_test (tc, test_init_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("zero");
  tcase_add_test (tc, test_zero_vec);
  suite_add_tcase (s, tc);

  return s;
}


#endif
