#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/Value.h>
#include <dolfin/function/Zero.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
struct E01 : public Value<E01>
{
  void eval(real * values, const real*) const
  {
    values[0] = 0.0;
  }
};

struct E12 : public Value<E12, 2>
{
  void eval(real * values, const real*) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
  }
};

struct E13 : public Value<E13, 3>
{
  void eval(real * values, const real*) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
    values[2] = 2.0;
  }
};

struct E22 : public Value<E22, 2, 2>
{
  void eval(real * values, const real*) const
  {
    values[0] = 0.0;
    values[1] = 1.0;
    values[2] = 2.0;
    values[3] = 3.0;
  }
};

struct E33 : public Value<E33, 3, 3>
{
  void eval(real * values, const real*) const
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
DOLFIN_START_TEST( test_Value )
  {
    {
      E01 e;
      e.disp();
    }
    //---
    {
      E12 e;
      e.disp();
    }
    //---
    {
      E13 e;
      e.disp();
    }
    //---
    {
      E22 e;
      e.disp();
    }
    //---
    {
      E33 e;
      e.disp();
    }
    //---
    {
      Zero<1> z1;
      z1.disp();
      Zero<2> z2;
      z2.disp();
      Zero<3> z3;
      z3.disp();
      Zero<2, 2> z22;
      z22.disp();
      Zero<3, 3> z33;
      z33.disp();
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
