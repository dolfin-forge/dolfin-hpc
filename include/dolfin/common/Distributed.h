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
  MPI::Communicator & comm();

  inline uint comm_rank() const;
  inline uint comm_size() const;

  inline bool distributed() const;

  /// Swap instances
  friend void swap( Distributed< T > & a, Distributed< T > & b )
  {
    using std::swap;
    swap( a.comm_, b.comm_ );
  }

protected:
  //----------------------------------------------------------------------------
  virtual ~Distributed();

  Distributed & operator=( Distributed const & other );

private:
  //----------------------------------------------------------------------------
  MPI::Communicator comm_;
};

//------------------------------------------------------------------------------
template < typename T >
Distributed< T >::Distributed( MPI::Communicator & comm )
  : comm_( DOLFIN_COMM_NULL )
{
#if HAVE_MPI
  /*
   * MPI 1.1:
   *  "A null handle argument is an erroneous IN argument in MPI calls"
   */
  if ( comm != DOLFIN_COMM_NULL )
    MPI::check_error( MPI_Comm_dup( comm, &comm_ ) );
#else
  MAYBE_UNUSED( comm );
#endif
}

//------------------------------------------------------------------------------
template < typename T >
Distributed< T >::Distributed( Distributed const & other )
  : comm_( DOLFIN_COMM_NULL )
{
  *this = other;
}
//------------------------------------------------------------------------------
template < typename T >
MPI::Communicator & Distributed< T >::comm()
{
  return comm_;
}

//------------------------------------------------------------------------------
template < typename T >
inline uint Distributed< T >::comm_rank() const
{
  int ret = 0;
#if HAVE_MPI
  if ( comm_ != DOLFIN_COMM_NULL )
    MPI::check_error( MPI_Comm_rank( comm_, &ret ) );
#endif
  return static_cast< uint >( ret );
}

//------------------------------------------------------------------------------
template < typename T >
inline uint Distributed< T >::comm_size() const
{
  int ret = 1;
#if HAVE_MPI
  if ( comm_ != DOLFIN_COMM_NULL )
    MPI::check_error( MPI_Comm_size( comm_, &ret ) );
#endif
  return static_cast< uint >( ret );
}

//------------------------------------------------------------------------------
template < typename T >
inline bool Distributed< T >::distributed() const
{
  return ( this->comm_size() > 1 );
}

//------------------------------------------------------------------------------
template < typename T >
Distributed< T >::~Distributed()
{
#if HAVE_MPI
  if ( comm_ != DOLFIN_COMM_NULL )
    MPI::check_error( MPI_Comm_free( &comm_ ) );
#endif
}

//------------------------------------------------------------------------------
template < typename T >
Distributed< T > & Distributed< T >::operator=( Distributed< T > const & other )
{
  if ( this != &other )
  {
#if HAVE_MPI
    if ( comm_ != DOLFIN_COMM_NULL )
      MPI::check_error( MPI_Comm_free( &comm_ ) );
    if ( other.comm_ != DOLFIN_COMM_NULL )
      MPI::check_error( MPI_Comm_dup( other.comm_, &comm_ ) );
#endif
  }
  return *this;
}

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_DISTRIBUTED_H_ */
