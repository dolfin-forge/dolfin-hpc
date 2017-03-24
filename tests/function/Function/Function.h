#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Test.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <check.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Function )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Function");
  {
    UnitInterval m0(16);
    ufl::FiniteElement DG0(ufl::Family::DG, m0.type(), 0);
    FiniteElementSpace Vh0(m0, DG0);
    Function U0(Vh0);

    //
    U0 = 1.0;
    ck_assert( U0.min() == 1.0 );
    ck_assert( U0.max() == 1.0 );
    U0.disp();


    UnitInterval m1(32);
    FiniteElementSpace Vh1(m1, DG0);
    Function U1(Vh1);

    //
    real * U1blck = U1.create_block();
    for (CellIterator cell(m1); !cell.end(); ++cell)
    {
      U1blck[cell->index()] = cell->index();
    }
    U1.set_block(U1blck);
    delete [] U1blck;

    //
    Function U2(m1);
    U2 = U1;

    U2.disp();



  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
