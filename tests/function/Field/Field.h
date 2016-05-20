#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Field.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Field )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Field");
  {
    UnitInterval m0(16);
    ufl::FiniteElement DG0(ufl::Family::DG, m0.type(), 0);
    FiniteElementSpace Vh0(m0, DG0);
    Field U0(Vh0);
  }
  T.end();
  //---
  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
