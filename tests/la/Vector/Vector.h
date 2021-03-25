
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
#ifdef HAVE_TRILINOS
#include <dolfin/la/trilinos/TrilinosVector.h>
DOLFIN_START_TEST( test_trilinos_vec )
{
  trilinos::Vector x1;
  x1.init(VEC_SIZE);
  trilinos::Vector x2( VEC_SIZE );
  trilinos::Vector x3( x2 );

  x1.disp();
  x1.zero();
  x1.size();
  x1.local_size();
  x1.offset();

  // global: get, set add
  // local: get, set add

  x1.axpy( 1.1, x2 );
  x1.norm( VectorNormType::l1 );
  x1.axpby( 1.1, x2, 2.3 );
  x1.norm( VectorNormType::l2 );
  x1.waxpy( 1.1, x2, x3 );
  x1.norm( VectorNormType::linf );
  x1.axpbypcz( 1.1, x2, 2.3, x3, 4.2 );

  x1.min();
  x1.max();

  x1 *= 3.2;
  x1 /= 4.2;
  x2  = 23.42;
  x1 *= x2;
  x1 += x2;
  x1 -= x2;
  x3  = x1;

  // x = 1.0;
  // x.zero();
  // ck_assert(x.max() == 0);

}DOLFIN_END_TEST
#else
DOLFIN_START_TEST( test_trilinos_vec )
{
}DOLFIN_END_TEST
#endif
//-----------------------------------------------------------------------------

#endif
