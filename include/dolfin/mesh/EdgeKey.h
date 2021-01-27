
#ifndef __DOLFIN_MESH_EDGE_KEY
#define __DOLFIN_MESH_EDGE_KEY

#include <dolfin/mesh/entities/Edge.h>

namespace dolfin
{

template < typename T = size_t >
struct EdgeKey : public std::pair< T, T >
{

  /// An edge contains a pair of vertices
  EdgeKey()
    : std::pair< T, T >()
  {
  }

  /// An edge contains a pair of vertices
  EdgeKey( Edge const & e )
    : std::pair< T, T >()
    , idx( std::rand() )
  {
    std::vector< size_t > const & v = e.entities( 0 );
    this->first  = static_cast< T >( v[0] < v[1] ? v[0] : v[1] );
    this->second = static_cast< T >( v[0] < v[1] ? v[1] : v[0] );
  }

  /// An edge contains a pair of vertices
  EdgeKey( T v0, T v1 )
    : std::pair< T, T >( v0 < v1 ? v0 : v1, v0 < v1 ? v1 : v0 )
    , idx( std::rand() )
  {
  }

  /// Construct a key from edge vertices
  inline void set( Edge const & e )
  {
    std::vector< size_t > const & v = e.entities( 0 );
    this->first  = static_cast< T >( v[0] < v[1] ? v[0] : v[1] );
    this->second = static_cast< T >( v[0] < v[1] ? v[1] : v[0] );
    idx          = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set( T const * v )
  {
    this->first  = ( v[0] < v[1] ? v[0] : v[1] );
    this->second = ( v[0] < v[1] ? v[1] : v[0] );
    idx          = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set( T v0, T v1 )
  {
    this->first  = ( v0 < v1 ? v0 : v1 );
    this->second = ( v0 < v1 ? v1 : v0 );
    idx          = std::rand();
  }

  ///
  inline auto hash() const -> size_t
  {
    return ( static_cast< std::size_t >( this->first )
             ^ static_cast< std::size_t >( this->second ) );
  }

  ///
  size_t idx { 0 };
};

} /* namespace dolfin */

namespace std
{

template < typename T >
struct hash< dolfin::EdgeKey< T > >
{
  inline auto operator()( dolfin::EdgeKey< T > const & e ) const -> std::size_t
  {
    return e.hash();
  }
};

} /* namespace std */

#endif /* __DOLFIN_MESH_EDGE_KEY */
