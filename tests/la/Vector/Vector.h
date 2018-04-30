#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/la/Vector.h>
#include <cstring>

#define VEC_SIZE 100

using namespace dolfin;
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_init_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  ck_assert(x.size() == VEC_SIZE);
  ck_assert(x.local_size() == VEC_SIZE);

}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_get_set_vec )
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
    ck_assert(data[i] == 3.1415);
  }

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_add_vec )
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
    ck_assert(data[i] == 2.0 * (3.1415));
  }

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_max_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = i;
  }

  x.set(data);

  ck_assert(x.max() == double (VEC_SIZE - 1));

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_min_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  double *data = new double[VEC_SIZE];

  for (int i = 0; i < VEC_SIZE; i++)
  {
    data[i] = i;
  }

  x.set(data);

  ck_assert(x.min() == 0.0);

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_op_assign_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    ck_assert(data[i] == 3.1457);
  }

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_op_mul_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x *= 2;

  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    ck_assert(data[i] == (2.0 * 3.1457));
  }

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_op_scale_vec )
{
  Vector x;
  x.init(VEC_SIZE);
  
  double *data = new double[VEC_SIZE];

  x = 3.1457;
  x /= 2;

  x.get(data);

  for (int i = 0; i < VEC_SIZE; i++)
  {
    ck_assert(data[i] == (3.1457 / 2.0));
  }

  delete [] data;
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_zero_vec )
{
  Vector x;
  x.init(VEC_SIZE);

  x = 1.0;
  x.zero();
  ck_assert(x.max() == 0);
  
}DOLFIN_END_TEST
//-----------------------------------------------------------------------------

DOLFIN_SUITE_BEGIN(test_suite_vec, "Vector")
{
  DOLFIN_TCASE_CREATE ("init");
  DOLFIN_TCASE_ADD(test_init_vec);

  DOLFIN_TCASE_CREATE ("get/set");
  DOLFIN_TCASE_ADD(test_get_set_vec);

  DOLFIN_TCASE_CREATE ("add");
  DOLFIN_TCASE_ADD(test_add_vec);

  DOLFIN_TCASE_CREATE ("max");
  DOLFIN_TCASE_ADD(test_max_vec);

  DOLFIN_TCASE_CREATE ("min");
  DOLFIN_TCASE_ADD(test_min_vec);

  DOLFIN_TCASE_CREATE ("operator =");
  DOLFIN_TCASE_ADD(test_op_assign_vec);

  DOLFIN_TCASE_CREATE ("operator *");
  DOLFIN_TCASE_ADD(test_op_mul_vec);

  DOLFIN_TCASE_CREATE ("operator /");
  DOLFIN_TCASE_ADD(test_op_scale_vec);

  DOLFIN_TCASE_CREATE ("zero");
  DOLFIN_TCASE_ADD(test_zero_vec);
}
DOLFIN_SUITE_END

#endif
