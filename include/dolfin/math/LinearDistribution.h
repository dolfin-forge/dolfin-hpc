
#ifndef __DOLFIN_LINEAR_DISTRIBUTION_H
#define __DOLFIN_LINEAR_DISTRIBUTION_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct LinearDistribution
{

  size_t global_size { 0 };
  size_t card { 0 };
  size_t rank { 0 };
  size_t L { 0 };
  size_t R { 0 };
  size_t offset { 0 };
  size_t size { 0 };

  ///
  LinearDistribution() = default;

  ///
  LinearDistribution( size_t global_size, size_t card, size_t rank )
    : global_size( 0 )
    , card( 0 )
    , rank( 0 )
    , L( 0 )
    , R( 0 )
    , offset( 0 )
    , size( 0 )
  {
    set( global_size, card, rank );
  }

  ///
  inline bool in_range( size_t index )
  {
    return ( ( offset <= index ) && ( index < offset + size ) );
  }

  ///
  inline size_t owner( size_t index ) const
  {
    return static_cast< size_t >( std::max(
      std::floor( ( real ) index / ( real )( L + 1 ) ),
      std::floor( ( real )( ( real ) index - ( real ) R ) / ( real ) L ) ) );
  }

  ///
  void set( size_t global_size, size_t card, size_t rank )
  {
    this->global_size = global_size, this->card = card;
    this->rank   = rank;
    this->L      = std::floor( ( real ) global_size / ( real ) card );
    this->R      = global_size % card;
    this->offset = rank * L + std::min( rank, R );
    this->size =
      std::floor( ( ( real ) global_size + ( real ) card - ( real ) rank - 1.0 )
                  / ( real ) card );
  }

  ///
  void disp() const
  {
    section( "LinearDistribution" );
    message( "global size    : %u", global_size );
    message( "num partitions : %u", card );
    message( "rank           : %u", rank );
    message( "offset         : %u", offset );
    message( "size           : %u", size );
    message( "quotient size  : %u", L );
    message( "remainder      : %u", R );
    end();
  }
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_LINEAR_DISTRIBUTION_H */
