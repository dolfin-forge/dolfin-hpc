#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/function/SpaceTimeFunction.h>
#include <dolfin/main/PE.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_SpaceTimeFunction )
{
  using namespace ElementLibrary;

  pushd( testresu( "function", "SpaceTimeFunction" ) );

  // create ufc::finite element
  ufc::finite_element * element = create_finite_element( elements[0] );

  // create ufc::dofmap
  ufc::dofmap * dofmap = create_dof_map( dofmaps[0] );

  Simplex< 1 > M;
  pushd( "0" );
  {
    SpaceTimeFunction S( "U", FiniteElementSpace( M, element, *dofmap ) );
    real const        dt = 0.1;
    for ( Time t( 0., 1. ); t.is_valid(); t += dt )
    {
      S( t ) = real( t );
      ck_assert( S.min() == S.max() );
      ck_assert( S.min() == real( t ) );
      S.save();
      // message("(%e, %e)", S.min(), S.max());
    }
    for ( Time t( 0., 1. ); t.is_valid(); t += dt )
    {
      S( t );
      S.eval();
      if ( PE::rank() == 0 )
      {
        message( "%e = (%e, %e)", real( t ), S.min(), S.max() );
      }
      ck_assert( S.min() == S.max() );
      ck_assert( S.min() == real( t ) );
    }
  }
  popd();
  pushd( "1" );
  {
    real const dt = 0.1;
    {
      SpaceTimeFunction S( "U", FiniteElementSpace( M, element, *dofmap ) );
      for ( Time t( 0., 1. ); t.is_valid(); t += dt )
      {
        S( t ) = real( t );
        ck_assert( S.min() == S.max() );
        ck_assert( S.min() == real( t ) );
        S.save();
        // message("(%e, %e)", S.min(), S.max());
      }
    }
    {
      SpaceTimeFunction S( "U", FiniteElementSpace( M, element, *dofmap ) );
      for ( Time t( 0., 1. ); t.is_valid(); t += dt )
      {
        S( t );
        S.eval();
        if ( PE::rank() == 0 )
        {
          message( "%e = (%e, %e)", real( t ), S.min(), S.max() );
        }
        ck_assert( S.min() == S.max() );
        ck_assert( S.min() == real( t ) );
      }
    }
  }
  popd();
  popd();

  delete element;
  delete dofmap;
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
