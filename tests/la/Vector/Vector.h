#ifdef HAVE_CHECK

#include <dolfin/la/Vector.h>
#include <cstring>

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
START_TEST( test_get_set_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = 3.1415;
  }

  x.set(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = 0.0;
  }
  
  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    fail_unless(data[i] == 3.1415);
  }

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_add_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = 3.1415;
  }

  x.set(data);
  x.add(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = 0.0;
  }

  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    fail_unless(data[i] == 2.0 * (3.1415));
  }

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_max_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = i;
  }

  x.set(data);

  fail_unless(x.max() == double (VEC_SIZE - 1));

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_min_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = i;
  }

  x.set(data);

  fail_unless(x.min() == 0.0);

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_op_assign_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    fail_unless(data[i] == 3.1457);
  }

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_op_mul_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x *= 2;

  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    fail_unless(data[i] == (2.0 * 3.1457));
  }

  delete [] data;
  
}END_TEST
//-----------------------------------------------------------------------------
START_TEST( test_op_scale_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x /= 2;

  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    fail_unless(data[i] == (3.1457 / 2.0));
  }

  delete [] data;
  
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

  tc = tcase_create ("get/set");
  tcase_add_test (tc, test_get_set_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("add");
  tcase_add_test (tc, test_add_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("max");
  tcase_add_test (tc, test_max_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("min");
  tcase_add_test (tc, test_min_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("operator =");
  tcase_add_test (tc, test_op_assign_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("operator *");
  tcase_add_test (tc, test_op_mul_vec);
  suite_add_tcase (s, tc);

  tc = tcase_create ("operator /");
  tcase_add_test (tc, test_op_scale_vec);
  suite_add_tcase (s, tc);

 tc = tcase_create ("zero");
  tcase_add_test (tc, test_zero_vec);
  suite_add_tcase (s, tc);

  return s;
}


#endif
