#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/log/log.h>
#include <dolfin/parameter/Parameter.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_parameter )
{
  message( "bool" );
  {
    parameter< bool > Q( false, Parameter::bool_t );
    message( "Q = %u", Q.get() );
    ck_assert( Q.get() == false );
    parameter< bool > R( true, Parameter::bool_t );
    message( "R = %u", R.get() );
    ck_assert( R.get() == true );
    ck_assert( Q.get() != R.get() );
  }
  //---
  message( "uint" );
  {
    parameter< uint > Q( 23, Parameter::uint_t );
    message( "Q = %u", Q.get() );
    ck_assert( Q.get() == 23 );
    parameter< uint > R( 42, Parameter::uint_t );
    message( "E = %u", R.get() );
    ck_assert( R.get() == 42 );
    ck_assert( Q.get() != R.get() );
  }
  //---
  message( "int" );
  {
    parameter< int > Q( 0, Parameter::int_t );
    message( "Q = %+d", Q.get() );
    ck_assert( Q.get() == 0 );
    parameter< int > R( 1, Parameter::int_t );
    message( "R = %+d", R.get() );
    ck_assert( R.get() == 1 );
    ck_assert( Q.get() != R.get() );
    parameter< int > S( -1, Parameter::int_t );
    message( "S = %+d", S.get() );
    ck_assert( S.get() == -1 );
    ck_assert( Q.get() != S.get() );
    ck_assert( R.get() != S.get() );
  }
  //---
  message( "real" );
  {
    parameter< real > Q( 0., Parameter::real_t );
    message( "Q = %+e", Q.get() );
    ck_assert( Q.get() == 0. );
    parameter< real > R( 1., Parameter::real_t );
    message( "R = %+e", R.get() );
    ck_assert( R.get() == 1. );
    ck_assert( Q.get() != R.get() );
    parameter< real > S( -1., Parameter::real_t );
    message( "S = %+e", S.get() );
    ck_assert( S.get() == -1. );
    ck_assert( Q.get() != S.get() );
    ck_assert( R.get() != S.get() );
  }
  //---
  message( "string" );
  {
    parameter< std::string > Q( "", Parameter::real_t );
    message( "Q = %s", Q.get().c_str() );
    ck_assert( Q.get() == "" );
    parameter< std::string > R( "1", Parameter::real_t );
    message( "R = %s", R.get().c_str() );
    ck_assert( R.get() == "1" );
    ck_assert( Q.get() != R.get() );
    parameter< std::string > S( "11", Parameter::real_t );
    message( "S = %s", S.get().c_str() );
    ck_assert( S.get() == "11" );
    ck_assert( Q.get() != S.get() );
    ck_assert( R.get() != S.get() );
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
