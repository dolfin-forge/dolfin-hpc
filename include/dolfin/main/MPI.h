// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MPI_H
#define __DOLFIN_MPI_H

#include <dolfin/common/assert.h>
#include <dolfin/common/types.h>
#include <dolfin/main/MPI_Datatypes.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>

#ifdef HAVE_MPI
#define DOLFIN_COMM_NULL MPI_COMM_NULL
#else
#define DOLFIN_COMM_NULL 0
#define MPI_SUCCESS 0
#endif

namespace dolfin
{

/// This class provides utility functions for easy access of the number of
/// processes and current process number.

class MPI
{

public:
#ifdef HAVE_MPI
  using Communicator = MPI_Comm;
  using offset_t     = MPI_Offset;
#else
  using Communicator = int;
  using offset_t     = long long;
#endif

  /*
   *  Local communicator
   */

  /// Return process rank in local communicator
  static auto rank() -> size_t;

  /// Return local communicator size
  static auto size() -> size_t;

  /// Return if the given rank is valid
  static auto is_valid_rank( size_t rank ) -> bool;

  /// Return if the current process is the root process in local communicator
  static auto is_root() -> bool;

  /// Return group identifier
  static auto group_id() -> size_t;

  /*
   *  Global communicator
   */

  /// Return process rank in global communicator
  static auto global_rank() -> size_t;

  /// Return global communicator size
  static auto global_size() -> size_t;

  /// Return if the current process is the root process in global communicator
  static auto is_global_root() -> bool;

  /// Return number of groups in global communicator
  static auto num_groups() -> size_t;

  /*
   *  Global communicator
   */

  /// Return seed value for current process
  static auto seed() -> size_t;

  ///
  static auto offset( size_t         xl,
                      size_t &       offset,
                      Communicator & comm = MPI::DOLFIN_COMM ) -> void;

  ////
  enum ReductionType
  {
    sum,
    min,
    max
  };

  //// Wrap in a template function to allow use of functors
  template < typename T >
  static auto bcast( T * x, int n, int r,
                     Communicator & comm = MPI::DOLFIN_COMM ) -> int;

  //// Wrap in a template function to allow use of functors
  template < typename T >
  static auto all_gather( T *            sendbuf,
                          int            sendcount,
                          T *            recvbuf,
                          int            recvcount,
                          Communicator & comm = MPI::DOLFIN_COMM ) -> int;

  //// Wrap in a template function to allow use of functors
  template < int Op, typename T >
  static auto all_reduce( T x, T & r, Communicator & comm = MPI::DOLFIN_COMM )
    -> int;

  //// Wrap in a template function to allow use of functors
  template < int Op, typename T >
  static auto all_reduce( T *            x,
                          T *            r,
                          size_t         count,
                          Communicator & comm = MPI::DOLFIN_COMM ) -> int;

  //// Wrap in a template function to allow use of functors
  template < int Op, typename T >
  static auto all_reduce_in_place( T &            r,
                                   Communicator & comm = MPI::DOLFIN_COMM )
    -> int;

  //// Wrap in a template function to allow use of functors
  template < int Op, typename T >
  static auto all_reduce_in_place( std::vector< T > & data,
                                   Communicator &     comm = MPI::DOLFIN_COMM )
    -> int;


  template < int Op, typename T >
  static auto reduce( T const *      sendbuf,
                      T *            recvbuf,
                      int            count,
                      int            root,
                      Communicator & comm = MPI::DOLFIN_COMM ) -> int;

  //// Wrap in a template function to allow use of functors
  template < typename T >
  static auto sendrecv( T *            sendbuf,
                        int            sendcount,
                        int            destination,
                        T *            recvbuf,
                        int            recvcount,
                        int            source,
                        int            tag,
                        Communicator & comm = MPI::DOLFIN_COMM ) -> int;

  template < typename T >
  static auto sendrecv( std::vector< T > & sendbuf,
                        int                destination,
                        std::vector< T > & recvbuf,
                        int                source,
                        int                tag,
                        Communicator &     comm = MPI::DOLFIN_COMM ) -> int;

  /// Start MPI timer
  static auto startTimer() -> void;

  /// Start MPI timer with external counter;
  static auto startTimer( dolfin::real & stime ) -> void;

  /// Stop MPI timer
  static auto stopTimer() -> real;

  /// Stop MPI timer
  static auto stopTimer( dolfin::real & stime ) -> real;

  /// Setup DOLFIN_COMM MPI communicator
  static auto initComm( int ngroups = 0 ) -> void;
  static auto finiComm() -> void;

  static Communicator DOLFIN_COMM_WORLD;
  static Communicator DOLFIN_COMM_SELF;
  static Communicator DOLFIN_COMM;

  /// Check for MPI errors
  static auto check_error( int const mpi_error ) -> int;

#ifdef HAVE_MPI

  static auto file_open( MPI_File &          file,
                         std::string const & filename,
                         int                 mode,
                         Communicator &      comm = MPI::DOLFIN_COMM,
                         MPI_Info            info = MPI_INFO_NULL ) -> void;

  template < typename T >
  static auto file_read_all( MPI_File &   file,
                             T &          element,
                             MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_read_all( MPI_File &   file,
                             T &          element,
                             size_t const size,
                             MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_read_at_all( MPI_File &   file,
                                T *          elements,
                                size_t const count,
                                offset_t     offset,
                                size_t const global_count = 0,
                                MPI_Status * status       = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_read_at_all( MPI_File &   file,
                                T &          element,
                                offset_t     offset,
                                MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_write_all( MPI_File &   file,
                              T const &    element,
                              MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_write_all( MPI_File &   file,
                              T const &    element,
                              size_t const size,
                              MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_write_at_all( MPI_File &   file,
                                 T const *    elements,
                                 size_t const count,
                                 offset_t     offset,
                                 size_t const global_count = 0,
                                 MPI_Status * status       = MPI_STATUS_IGNORE )
    -> offset_t;

  template < typename T >
  static auto file_write_at_all( MPI_File &   file,
                                 T const &    element,
                                 offset_t     offset,
                                 MPI_Status * status = MPI_STATUS_IGNORE )
    -> offset_t;

  static auto file_close( MPI_File & file ) -> void;

  static auto check_file_status( int const mpi_error ) -> void;

#endif

  template < typename T >
  static auto exscan_sum( T const *      send,
                          T *            recv,
                          int            count,
                          Communicator & comm = MPI::DOLFIN_COMM ) -> void;

private:
  static real time_;
  static bool init_;

  using Context = struct
  {
    int global_rank;
    int global_size;
    int group_cnt;
    int group_idx;
    int rank;
    int size;
    int seed;
  };

  static Context ctx_;
};

// Defines default type to MPI communication domain.

using Comm = MPI::Communicator;

#define DOLFIN_COMM_WORLD MPI::DOLFIN_COMM_WORLD
#define DOLFIN_COMM_SELF MPI::DOLFIN_COMM_SELF
#define DOLFIN_COMM MPI::DOLFIN_COMM

//-----------------------------------------------------------------------------

#if HAVE_MPI

//-----------------------------------------------------------------------------
template < typename T >
inline auto MPI::bcast( T * x, int n, int r, Communicator & comm ) -> int
{
  return MPI::check_error( MPI_Bcast( x, n, MPI_type< T >::value, r, comm ) );
}
//-----------------------------------------------------------------------------
template < typename T >
inline auto MPI::all_gather( T *            sendbuf,
                             int            sendcount,
                             T *            recvbuf,
                             int            recvcount,
                             Communicator & comm ) -> int
{
  return MPI::check_error( MPI_Allgather( sendbuf,
                                          sendcount,
                                          MPI_type< T >::value,
                                          recvbuf,
                                          recvcount,
                                          MPI_type< T >::value,
                                          comm ) );
}
//-----------------------------------------------------------------------------
// unfortunately c++ does not allow partial function template specialization
// so we have to move the  template to a mpi_helper class
template < int Op, typename T >
struct mpi_helper
{
  inline static auto reduce( T const * sendbuf, T * recvbuf, int count,
                             int root, MPI::Communicator & comm ) -> int;

  inline static auto all_reduce( T * x, T * r, int count,
                                 MPI::Communicator & comm ) -> int;

  inline static auto all_reduce_in_place( T * r, int count,
                                          MPI::Communicator & comm ) -> int;
};

template < typename T >
struct mpi_helper< MPI::sum, T >
{
  inline static auto reduce( T const * sendbuf, T * recvbuf, int count,
                             int root, MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Reduce( sendbuf, recvbuf, count,
                                         MPI_type< T >::value,
                                         MPI_SUM, root, comm ) );
  }

  inline static auto all_reduce( T * x, T * r, int count,
                                 MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( x, r, count,
                                            MPI_type< T >::value,
                                            MPI_SUM, comm ) );
  }

  inline static auto all_reduce_in_place( T * r, int count,
                                          MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( MPI_IN_PLACE, r, count,
                                            MPI_type< T >::value,
                                            MPI_SUM, comm ) );
  }
};

template < typename T >
struct mpi_helper< MPI::min, T >
{
  inline static auto reduce( T const * sendbuf, T * recvbuf, int count,
                             int root, MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Reduce( sendbuf, recvbuf, count,
                                         MPI_type< T >::value,
                                         MPI_MIN, root, comm ) );
  }

  inline static auto all_reduce( T * x, T * r, int count,
                                 MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( x, r, count,
                                            MPI_type< T >::value,
                                            MPI_MIN, comm ) );
  }

  inline static auto all_reduce_in_place( T * r, int count,
                                          MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( MPI_IN_PLACE, r, count,
                                            MPI_type< T >::value,
                                            MPI_MIN, comm ) );
  }
};

template < typename T >
struct mpi_helper< MPI::max, T >
{
  inline static auto reduce( T const * sendbuf, T * recvbuf, int count,
                             int root, MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Reduce( sendbuf, recvbuf, count,
                                         MPI_type< T >::value,
                                         MPI_MAX, root, comm ) );
  }

  inline static auto all_reduce( T * x, T * r, int count,
                                 MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( x, r, count,
                                            MPI_type< T >::value,
                                            MPI_MAX, comm ) );
  }

  inline static auto all_reduce_in_place( T * r, int count,
                                          MPI::Communicator & comm ) -> int
  {
    return MPI::check_error( MPI_Allreduce( MPI_IN_PLACE, r, count,
                                            MPI_type< T >::value,
                                            MPI_MAX, comm ) );
  }
};


template < int Op, typename T >
inline auto MPI::reduce( T const * sendbuf, T * recvbuf,
                         int count, int root, Communicator & comm ) -> int
{
  return mpi_helper< Op, T >::reduce( sendbuf, recvbuf, count, root, comm );
}

template < int Op, typename T >
inline auto MPI::all_reduce( T x, T & r, Communicator & comm ) -> int
{
  return mpi_helper< Op, T >::all_reduce( &x, &r, 1, comm );
}

template < int Op, typename T >
inline auto MPI::all_reduce( T * x, T * r, size_t count, Communicator & comm )
  -> int
{
  return mpi_helper< Op, T >::all_reduce( x, r, count, comm );
}

template < int Op, typename T >
inline auto MPI::all_reduce_in_place( T & r, Communicator & comm ) -> int
{
  return mpi_helper< Op, T >::all_reduce_in_place( &r, 1, comm );
}

template < int Op, typename T >
inline auto MPI::all_reduce_in_place( std::vector< T > & data,
                                      Communicator &     comm ) -> int
{
  return mpi_helper< Op, T >::all_reduce_in_place( data.data(), data.size(), comm );
}

//-----------------------------------------------------------------------------
template <>
inline auto MPI::sendrecv( bool *         sendbuf,
                           int            sendcount,
                           int            destination,
                           bool *         recvbuf,
                           int            recvcount,
                           int            source,
                           int            tag,
                           Communicator & comm ) -> int
{
  MPI_Status status;
  int        recv_count;
  int        send_bool = sendcount * sizeof( bool );
  int        recv_bool = recvcount * sizeof( bool );
  MPI::check_error( MPI_Sendrecv( static_cast< void * >( sendbuf ),
                                  send_bool,
                                  MPI_BYTE,
                                  destination,
                                  tag,
                                  static_cast< void * >( recvbuf ),
                                  recv_bool,
                                  MPI_BYTE,
                                  source,
                                  tag,
                                  comm,
                                  &status ) );
  MPI::check_error( MPI_Get_count( &status, MPI_BYTE, &recv_count ) );
  return recv_count / sizeof( bool );
}

template < typename T >
inline auto MPI::sendrecv( T *            sendbuf,
                           int            sendcount,
                           int            destination,
                           T *            recvbuf,
                           int            recvcount,
                           int            source,
                           int            tag,
                           Communicator & comm ) -> int
{
  MPI_Status status;
  int        recv_count;
  MPI::check_error( MPI_Sendrecv( static_cast< void * >( sendbuf ),
                                  sendcount,
                                  MPI_type< T >::value,
                                  destination,
                                  tag,
                                  static_cast< void * >( recvbuf ),
                                  recvcount,
                                  MPI_type< T >::value,
                                  source,
                                  tag,
                                  comm,
                                  &status ) );
  MPI::check_error(
    MPI_Get_count( &status, MPI_type< T >::value, &recv_count ) );
  return recv_count;
}

template < typename T >
inline auto MPI::sendrecv( std::vector< T > & sendbuf,
                           int                destination,
                           std::vector< T > & recvbuf,
                           int                source,
                           int                tag,
                           Communicator &     comm ) -> int
{
  return MPI::sendrecv< T >( sendbuf.data(),
                             sendbuf.size(),
                             destination,
                             recvbuf.data(),
                             recvbuf.size(),
                             source,
                             tag,
                             comm );
}
//-----------------------------------------------------------------------------
#if not defined( DEBUG )
inline int MPI::check_error( int const mpi_error )
{
  return mpi_error;
}
#endif
//-----------------------------------------------------------------------------

template < typename T >
auto MPI::file_read_all( MPI_File & file, T & element, MPI_Status * status )
  -> MPI::offset_t
{
  check_error( MPI_File_read_all( file,
                                  static_cast< void * >( &element ),
                                  1,
                                  MPI_type< T >::value,
                                  status ) );
  return static_cast< offset_t >( sizeof( T ) );
}

template < typename T >
auto MPI::file_read_all( MPI_File &   file,
                         T &          element,
                         size_t const size,
                         MPI_Status * status ) -> MPI::offset_t
{
  check_error( MPI_File_read_all(
    file, static_cast< void * >( &element ), size, MPI_BYTE, status ) );
  return static_cast< offset_t >( size );
}

//-----------------------------------------------------------------------------

template < typename T >
auto MPI::file_read_at_all( MPI_File &   file,
                            T *          elements,
                            size_t const count,
                            offset_t     offset,
                            size_t const global_count,
                            MPI_Status * status ) -> MPI::offset_t
{
  dolfin_assert( count <= global_count or global_count == 0 );
  check_error( MPI_File_read_at_all( file,
                                     offset,
                                     static_cast< void * >( elements ),
                                     count,
                                     MPI_type< T >::value,
                                     status ) );
  return static_cast< offset_t >( global_count * sizeof( T ) );
}

template < typename T >
auto MPI::file_read_at_all( MPI_File &   file,
                            T &          element,
                            offset_t     offset,
                            MPI_Status * status ) -> MPI::offset_t
{
  check_error( MPI_File_read_at_all( file,
                                     offset,
                                     static_cast< void * >( &element ),
                                     sizeof( T ),
                                     MPI_BYTE,
                                     status ) );
  return static_cast< offset_t >( sizeof( T ) );
}

//-----------------------------------------------------------------------------

template < typename T >
auto MPI::file_write_all( MPI_File &   file,
                          T const &    element,
                          MPI_Status * status ) -> MPI::offset_t
{
  check_error( MPI_File_write_all( file,
                                   static_cast< void const * >( &element ),
                                   1,
                                   MPI_type< T >::value,
                                   status ) );
  return static_cast< offset_t >( sizeof( T ) );
}

template < typename T >
auto MPI::file_write_all( MPI_File &   file,
                          T const &    element,
                          size_t const size,
                          MPI_Status * status ) -> MPI::offset_t
{
  check_error( MPI_File_write_all(
    file, static_cast< void const * >( &element ), size, MPI_BYTE, status ) );
  return static_cast< offset_t >( size );
}

//-----------------------------------------------------------------------------

template < typename T >
auto MPI::file_write_at_all( MPI_File &   file,
                             T const *    elements,
                             size_t const count,
                             offset_t     offset,
                             size_t const global_count,
                             MPI_Status * status ) -> MPI::offset_t
{
  dolfin_assert( count <= global_count or global_count == 0 );
  check_error( MPI_File_write_at_all( file,
                                      offset,
                                      static_cast< void const * >( elements ),
                                      count,
                                      MPI_type< T >::value,
                                      status ) );
  return static_cast< offset_t >( global_count * sizeof( T ) );
}

template < typename T >
auto MPI::file_write_at_all( MPI_File &   file,
                             T const &    element,
                             offset_t     offset,
                             MPI_Status * status ) -> MPI::offset_t
{
  check_error( MPI_File_write_at_all( file,
                                      offset,
                                      static_cast< void const * >( &element ),
                                      sizeof( T ),
                                      MPI_BYTE,
                                      status ) );
  return static_cast< offset_t >( sizeof( T ) );
}

//-----------------------------------------------------------------------------

template < typename T >
void MPI::exscan_sum( T const * send, T * recv, int count, Communicator & comm )
{
#if ( MPI_VERSION > 1 )
  MPI::check_error( MPI_Exscan( static_cast< void const * >( send ),
                                static_cast< void * >( recv ),
                                count,
                                MPI_type< T >::value,
                                MPI_SUM,
                                comm ) );
#else
  MPI::check_error( MPI_Scan( static_cast< void const * >( send ),
                              static_cast< void * >( recv ),
                              count,
                              MPI_type< T >::value,
                              MPI_SUM,
                              comm ) );

  for ( int i = 0; i < count; ++i )
    recv[i] -= send[i];
#endif
}

#else

//-----------------------------------------------------------------------------
template < typename T >
inline int MPI::bcast( T *, int, int, Communicator & )
{
  return MPI_SUCCESS;
}

//-----------------------------------------------------------------------------
template < typename T >
inline int MPI::all_gather( T *, int, T *, int, Communicator & )
{
  return MPI_SUCCESS;
}

//-----------------------------------------------------------------------------
template < int Op, typename T >
inline int MPI::all_reduce( T x, T & r, Communicator & )
{
  r = x;
  return MPI_SUCCESS;
}

//-----------------------------------------------------------------------------
template < typename T >
inline int MPI::sendrecv( T *, int, int, T *, int, int, int, Communicator & )
{
  return MPI_SUCCESS;
}

//-----------------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_H */
