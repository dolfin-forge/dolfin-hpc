#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/common/Array.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Connectivity.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template < class T, size_t D >
void check_regular_connectivity()
{
  message( "Connectivity: check %u-regular connectivity", D );
  for ( size_t i = 0; i <= 4; ++i )
  {
    // Create connectivity and fill
    Connectivity C( 1 << i, D );
    ck_assert_uint_eq( C.entries(), ( 1 << i ) * D );

    size_t ei = 0;
    for ( size_t it = 0; it < C.order(); ++it )
    {
      if ( D > 0 )
        ck_assert( not C[it].empty() );
      else
        ck_assert( C[it].empty() );

      for ( size_t cit = 0; cit < C[it].size(); ++cit )
        C[it][cit] = ei++;
    }
    // Check entries
    for ( size_t vi = 0; vi < C.order(); ++vi )
    {
      ei = 0;
      for ( size_t cvi = 0; cvi < C.degree( vi ); ++cvi, ++ei )
      {
        ck_assert( C[vi][cvi] == ( vi * D + ei ) );
      }
    }

    // Perform left remapping
    {
      Array< size_t > L( C.order() );
      if ( C.order() )
        range( L.data(), L.data() + L.size(), L.size() - 1, -1 );
      C.remap_l( L );
      size_t vj = C.order();
      for ( size_t vi = 0; vi < C.order(); ++vi )
      {
        ei = 0;
        --vj;
        for ( size_t cvi = 0; cvi < C.degree( vi ); ++cvi, ++ei )
        {
          ck_assert( C[vi][cvi] == ( vj * D + ei ) );
        }
      }
    }

    // Perform right remapping
    {
      Array< size_t > R( C.entries() );
      if ( C.entries() )
        range( R.data(), R.data() + R.size(), R.size() - 1, -1 );
      C.remap_r( R );
      for ( size_t vi = 0; vi < C.order(); ++vi )
      {
        ei = 1;
        for ( size_t cvi = 0; cvi < C.degree( vi ); ++cvi, ++ei )
        {
          ck_assert( C[vi][cvi] == ( vi * D + D - ei ) );
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Connectivity )
{
  //--- Create empty connectivity, copy, and iterate
  {
    Connectivity C0( 0, 0 );
    ck_assert( C0.order() == 0 );
    ck_assert( C0.min_degree() == 0 );
    ck_assert( C0.max_degree() == 0 );
    ck_assert( C0.entries() == 0 );
    ck_assert( C0().empty() == true );
    Connectivity C1( C0 );
    ck_assert( C1().empty() == true );
    ck_assert( C0 == C1 );
    for ( size_t it = 0; it < C1.order(); ++it )
    {
      error( "Iteration on empty connectivity" );
    }
    // Empty remappings
    Array< size_t > L;
    C0.remap_l( L );
    Array< size_t > R;
    C0.remap_r( R );
  }
  //--- Create 2-connectivity, verify basic data, and remap
  {
    Connectivity C( 4, 2 );
    ck_assert_uint_eq( C.order(), 4 );
    ck_assert_uint_eq( C.min_degree(), 2 );
    ck_assert_uint_eq( C.max_degree(), 2 );
    ck_assert_uint_eq( C.entries(), 8 );

    C[0][0] = 0;
    C[0][1] = 1;
    C[1][0] = 1;
    C[1][1] = 2;
    C[2][0] = 2;
    C[2][1] = 3;
    C[3][0] = 3;
    C[3][1] = 4;

    ck_assert_uint_eq( C[0][0], 0 );
    ck_assert_uint_eq( C[0][1], 1 );
    ck_assert_uint_eq( C[1][0], 1 );
    ck_assert_uint_eq( C[1][1], 2 );
    ck_assert_uint_eq( C[2][0], 2 );
    ck_assert_uint_eq( C[2][1], 3 );
    ck_assert_uint_eq( C[3][0], 3 );
    ck_assert_uint_eq( C[3][1], 4 );

    Array< Array< size_t > > connect;
    connect << C;
    Connectivity D( connect );
    ck_assert( C == D );
    Array< size_t > L( 4 );
    L[0] = 3;
    L[1] = 2;
    L[2] = 1;
    L[3] = 0;
    C.remap_l( L );

    ck_assert_uint_eq( C[3][0], 0 );
    ck_assert_uint_eq( C[3][1], 1 );
    ck_assert_uint_eq( C[2][0], 1 );
    ck_assert_uint_eq( C[2][1], 2 );
    ck_assert_uint_eq( C[1][0], 2 );
    ck_assert_uint_eq( C[1][1], 3 );
    ck_assert_uint_eq( C[0][0], 3 );
    ck_assert_uint_eq( C[0][1], 4 );

    Array< size_t > R( 5 );
    R[0] = 4;
    R[1] = 3;
    R[2] = 2;
    R[3] = 1;
    R[4] = 0;
    C.remap_r( R );

    ck_assert_uint_eq( C[3][0], 4 );
    ck_assert_uint_eq( C[3][1], 3 );
    ck_assert_uint_eq( C[2][0], 3 );
    ck_assert_uint_eq( C[2][1], 2 );
    ck_assert_uint_eq( C[1][0], 2 );
    ck_assert_uint_eq( C[1][1], 1 );
    ck_assert_uint_eq( C[0][0], 1 );
    ck_assert_uint_eq( C[0][1], 0 );
  }
  //--- Create k-connectivities and perform basic tests
  {
    check_regular_connectivity< Connectivity, 0 >();
    check_regular_connectivity< Connectivity, 1 >();
    check_regular_connectivity< Connectivity, 2 >();
    check_regular_connectivity< Connectivity, 3 >();
    check_regular_connectivity< Connectivity, 4 >();
    check_regular_connectivity< Connectivity, 5 >();
    check_regular_connectivity< Connectivity, 6 >();
    check_regular_connectivity< Connectivity, 8 >();
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
