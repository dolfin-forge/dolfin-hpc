#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Value.h>

using namespace dolfin;

struct E01 : public Value<E01>
{
  void eval(real * values, const real* x) const
  {
    values[0] = 0.0;
  }
};

struct E12 : public Value<E12,2>
{
  void eval(real * values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
  }
};

struct E13 : public Value<E13,3>
{
  void eval(real * values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
    values[2] = 2.0;
  }
};

struct E22 : public Value<E22,2, 2>
{
  void eval(real * values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
    values[2] = 2.0;
    values[3] = 3.0;
  }
};

struct E33 : public Value<E33,3, 3>
{
  void eval(real * values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
    values[2] = 2.0;
    values[3] = 3.0;
    values[4] = 4.0;
    values[5] = 5.0;
    values[3] = 6.0;
    values[4] = 7.0;
    values[5] = 8.0;
  }
};

//-----------------------------------------------------------------------------
START_TEST( test_Value )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Value : E01");
  {
    E01 e;
    e.disp();
  }
  T.end();
  //---
  T.begin("test_Value : E12");
  {
    E12 e;
    e.disp();
  }
  T.end();
  //---
  T.begin("test_Value : E13");
  {
    E13 e;
    e.disp();
  }
  T.end();
  //---
  T.begin("test_Value : E22");
  {
    E22 e;
    e.disp();
  }
  T.end();
  //---
  T.begin("test_Value : E33");
  {
    E33 e;
    e.disp();
  }
  T.end();
  //---
  T.begin("test_Value : Zero");
  {
    Zero<1> z1; z1.disp();
    Zero<2> z2; z2.disp();
    Zero<3> z3; z3.disp();
    Zero<2, 2> z22; z22.disp();
    Zero<3, 3> z33; z33.disp();
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
