#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/unitmeshes/UnitInterval.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Function )
{
  using namespace ElementLibrary;
  {
    // create ufc::finite element
    ufc::finite_element * element = create_finite_element( elements[0] );

    // create ufc::dofmap
    ufc::dofmap * dofmap = create_dof_map( dofmaps[0] );

    UnitInterval       m0( 16 );
    FiniteElementSpace Vh0( m0, element, *dofmap );
    Function           U0( Vh0 );

    //---
    U0 = 1.0;
    ck_assert( U0.min() == 1.0 );
    ck_assert( U0.max() == 1.0 );

    UnitInterval       m1( 32 );
    FiniteElementSpace Vh1( m1, element, *dofmap );
    Function           U1( Vh1 );

    //---
    real * U1blck = U1.create_block();
    for ( CellIterator cell( m1 ); !cell.end(); ++cell )
    {
      U1blck[cell->index()] = cell->index();
    }
    U1.set_block( U1blck );
    delete[] U1blck;

    //---
    Function U2( m1 );
    U2 = U1;

    delete element;
    delete dofmap;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
