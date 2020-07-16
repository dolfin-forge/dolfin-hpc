
#include <dolfin/common/AdjacentMapping.h>

#include <dolfin/common/SharedIterator.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SharedMapping::SharedMapping( DistributedData const & data )
  : data_( data )
  , mappings_()
  , send_min_( 0 )
  , send_max_( 0 )
{
  if ( !data_.is_finalized() )
  {
    error( "SharedMapping : distributed data is not finalized" );
  }

#if HAVE_MPI

  // Collect entities by adjacent rank
  for ( SharedIterator it( data ); it.valid(); ++it )
  {
    for( uint const & id : it.adj() )
      mappings_[id].send.push_back( it.global_index() );
  }

  dolfin_assert( mappings_.size() == data.get_adj_ranks().size() );
  send_max_ = 0;
  send_min_ = data.num_shared();

  //
  Array< MPI_Request> sendreq( mappings_.size() );
  Array< MPI_Request> recvreq( mappings_.size() );
  Array< MPI_Status>  status( mappings_.size() );

  uint i = 0;
  for ( std::pair< uint const, AdjacentMapping > & it : mappings_ )
  {
    dolfin_assert( it.first != MPI::rank() );

    AdjacentMapping & amap = it.second;

    // Update bounds
    send_max_ = std::max( send_max_, ( uint ) amap.send.size() );
    send_min_ = std::min( send_min_, ( uint ) amap.send.size() );
    //
    MPI::check_error( MPI_Isend( amap.send.data(), amap.send.size(),
                                 MPI_UNSIGNED, it.first, 0,
                                 data.comm(), &sendreq[i] ) );
    // Resize buffer
    amap.recv.resize( amap.send.size() );
    MPI::check_error( MPI_Irecv( amap.recv.data(), amap.recv.size(),
                                 MPI_UNSIGNED, it.first, 0,
                                 data.comm(), &recvreq[i] ) );

    ++i;
  }

  for ( std::pair< uint const, AdjacentMapping > & it : mappings_ )
  {
    AdjacentMapping & amap = it.second;
    data_.get_local( amap.send.size(), amap.send.data(), amap.send.data() );
  }

  MPI::check_error( MPI_Waitall( mappings_.size(), &sendreq[0], &status[0] ) );

  i = 0;
  for ( std::pair< uint const, AdjacentMapping > & it : mappings_ )
  {
    AdjacentMapping & amap = it.second;
    int recvcount = 0;

    MPI::check_error( MPI_Wait( &recvreq[i], &status[i] ) );
    MPI::check_error( MPI_Get_count( &status[i], MPI_UNSIGNED, &recvcount ) );

    if ( static_cast< uint >( recvcount ) != amap.recv.size() )
    {
      error( "AdjacentMapping : inconsistent count %u from rank %u: expected %u",
             MPI::rank(), recvcount, amap.recv.size() );
    }

    data_.get_local( amap.recv.size(), amap.recv.data(), amap.recv.data() );
    ++i;
  }

#endif /* HAVE_MPI */
}
//-----------------------------------------------------------------------------
SharedMapping::SharedMapping( SharedMapping const & other )
  : data_( other.data_ )
  , mappings_( other.mappings_ )
  , send_min_( other.send_min_ )
  , send_max_( other.send_max_ )
{
}
//-----------------------------------------------------------------------------
SharedMapping::~SharedMapping()
{
}
//-----------------------------------------------------------------------------
Array< uint > const & SharedMapping::to( uint rank ) const
{
  _map< uint, AdjacentMapping >::const_iterator it = mappings_.find( rank );
  if ( it == mappings_.end() )
  {
    error( "SharedMapping : invalid adjacent %u", rank );
  }
  return it->second.send;
}
//-----------------------------------------------------------------------------
Array< uint > const & SharedMapping::from( uint rank ) const
{
  _map< uint, AdjacentMapping >::const_iterator it = mappings_.find( rank );
  if ( it == mappings_.end() )
  {
    error( "SharedMapping : invalid adjacent %u", rank );
  }
  return it->second.recv;
}
//-----------------------------------------------------------------------------
void SharedMapping::disp() const
{
  section( "SharedMapping" );
  message( "number of adjacents : %u", mappings_.size() );
  message( "minimum size        : %u", send_min_ );
  message( "maximum size        : %u", send_max_ );
  end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
