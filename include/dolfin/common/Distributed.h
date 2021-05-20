// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_DISTRIBUTED_H_
#define __DOLFIN_COMMON_DISTRIBUTED_H_

#include <dolfin/common/maybe_unused.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

template < class T >
class Distributed
{
public:
  //----------------------------------------------------------------------------

  Distributed( MPI::Communicator & comm );

  Distributed( Distributed const & other );

  //----------------------------------------------------------------------------

  // access data
  auto comm() -> MPI::Communicator &;
  auto comm() const -> MPI::Communicator const &;

  inline auto comm_rank() const -> size_t;
  inline auto comm_size() const -> size_t;

  inline auto distributed() const -> bool;

  /// Swap instances
  friend void swap( Distributed< T > & a, Distributed< T > & b )
  {
    using std::swap;
    swap( a.comm_, b.comm_ );
  }

protected:
  //----------------------------------------------------------------------------
  virtual ~Distributed();

  auto operator=( Distributed const & other ) -> Distributed &;

private:
  //----------------------------------------------------------------------------
  MPI::Communicator comm_;
  size_t            comm_rank_;
  size_t            comm_size_;
};

//------------------------------------------------------------------------------

template < typename T >
Distributed< T >::Distributed( MPI::Communicator & comm )
  : comm_( DOLFIN_COMM_NULL )
  , comm_rank_( 0 )
  , comm_size_( 1 )
{
#if DOLFIN_HAVE_MPI
  /*
   * MPI 1.1:
   *  "A null handle argument is an erroneous IN argument in MPI calls"
   */
  if ( comm != DOLFIN_COMM_NULL )
  {
    MPI::check_error( MPI_Comm_dup( comm, &comm_ ) );

    int ret = 0;
    MPI::check_error( MPI_Comm_rank( comm_, &ret ) );
    comm_rank_ = static_cast< size_t >( ret );

    ret = 0;
    MPI::check_error( MPI_Comm_size( comm_, &ret ) );
    comm_size_ = static_cast< size_t >( ret );
  }
#else
  MAYBE_UNUSED( comm );
#endif
}

//------------------------------------------------------------------------------

template < typename T >
Distributed< T >::Distributed( Distributed const & other )
  : comm_( DOLFIN_COMM_NULL )
  , comm_rank_( 0 )
  , comm_size_( 1 )
{
  *this = other;
}
//------------------------------------------------------------------------------

template < typename T >
auto Distributed< T >::comm() -> MPI::Communicator &
{
  return comm_;
}
//------------------------------------------------------------------------------

template < typename T >
auto Distributed< T >::comm() const -> MPI::Communicator const &
{
  return comm_;
}

//------------------------------------------------------------------------------

template < typename T >
inline auto Distributed< T >::comm_rank() const -> size_t
{
  return comm_rank_;
}

//------------------------------------------------------------------------------

template < typename T >
inline auto Distributed< T >::comm_size() const -> size_t
{
  return comm_size_;
}

//------------------------------------------------------------------------------

template < typename T >
inline auto Distributed< T >::distributed() const -> bool
{
  return comm_size_ > 1;
}

//------------------------------------------------------------------------------

template < typename T >
Distributed< T >::~Distributed()
{
#if DOLFIN_HAVE_MPI
  if ( comm_ != DOLFIN_COMM_NULL )
    MPI::check_error( MPI_Comm_free( &comm_ ) );
#endif
}

//------------------------------------------------------------------------------

template < typename T >
auto Distributed< T >::operator=( Distributed< T > const & other )
  -> Distributed< T > &
{
  if ( this != &other )
  {
#if DOLFIN_HAVE_MPI
    if ( comm_ != DOLFIN_COMM_NULL )
      MPI::check_error( MPI_Comm_free( &comm_ ) );
    if ( other.comm_ != DOLFIN_COMM_NULL )
    {
      MPI::check_error( MPI_Comm_dup( other.comm_, &comm_ ) );

      int ret = 0;
      MPI::check_error( MPI_Comm_rank( comm_, &ret ) );
      comm_rank_ = static_cast< size_t >( ret );

      ret = 0;
      MPI::check_error( MPI_Comm_size( comm_, &ret ) );
      comm_size_ = static_cast< size_t >( ret );
    }
#endif
  }
  return *this;
}

//------------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_DISTRIBUTED_H_ */
