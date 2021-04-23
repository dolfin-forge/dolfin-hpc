
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

#if DOLFIN_HAVE_MPI

  // Collect entities by adjacent rank
  for ( SharedIterator it( data ); it.valid(); ++it )
  {
    for ( size_t const & id : it.adj() )
      mappings_[id].send.push_back( it.global_index() );
  }

  dolfin_assert( mappings_.size() == data.get_adj_ranks().size() );
  send_max_ = 0;
  send_min_ = data.num_shared();

  //
  std::vector< MPI_Request > sendreq( mappings_.size() );
  std::vector< MPI_Request > recvreq( mappings_.size() );
  std::vector< MPI_Status >  status( mappings_.size() );

  size_t i = 0;
  for ( std::pair< size_t const, AdjacentMapping > & it : mappings_ )
  {
    dolfin_assert( it.first != MPI::rank() );

    AdjacentMapping & amap = it.second;

    // Update bounds
    send_max_ = std::max( send_max_, amap.send.size() );
    send_min_ = std::min( send_min_, amap.send.size() );
    //
    MPI::check_error( MPI_Isend( amap.send.data(), amap.send.size(),
                                 MPI_type< size_t >::value,
                                 it.first, 0, data.comm(), &sendreq[i] ) );
    // Resize buffer
    amap.recv.resize( amap.send.size() );
    MPI::check_error( MPI_Irecv( amap.recv.data(), amap.recv.size(),
                                 MPI_type< size_t >::value,
                                 it.first, 0, data.comm(), &recvreq[i] ) );

    ++i;
  }

  for ( std::pair< size_t const, AdjacentMapping > & it : mappings_ )
  {
    AdjacentMapping & amap = it.second;
    data_.get_local( amap.send.size(), amap.send.data(), amap.send.data() );
  }

  MPI::check_error( MPI_Waitall( mappings_.size(), &sendreq[0], &status[0] ) );

  i = 0;
  for ( std::pair< size_t const, AdjacentMapping > & it : mappings_ )
  {
    AdjacentMapping & amap      = it.second;
    int               recvcount = 0;

    MPI::check_error( MPI_Wait( &recvreq[i], &status[i] ) );
    MPI::check_error( MPI_Get_count( &status[i], MPI_type< size_t >::value, &recvcount ) );

    if ( static_cast< size_t >( recvcount ) != amap.recv.size() )
    {
      error( "AdjacentMapping : inconsistent count %u from rank %u: expected %u",
             MPI::rank(), recvcount, amap.recv.size() );
    }

    data_.get_local( amap.recv.size(), amap.recv.data(), amap.recv.data() );
    ++i;
  }

  // init empty AdjacentMappings if they dont exist
  for ( size_t rank = 0; rank < MPI::size(); ++rank )
    if ( mappings_.find( rank ) == mappings_.end() )
      mappings_[rank] = AdjacentMapping();

#endif /* DOLFIN_HAVE_MPI */
}

//-----------------------------------------------------------------------------

auto SharedMapping::to( size_t rank ) const -> std::vector< size_t > const &
{
  _map< size_t, AdjacentMapping >::const_iterator it = mappings_.find( rank );
  if ( it == mappings_.end() )
  {
    error( "SharedMapping : invalid adjacent %u", rank );
  }
  return it->second.send;
}

//-----------------------------------------------------------------------------

auto SharedMapping::from( size_t rank ) const -> std::vector< size_t > const &
{
  _map< size_t, AdjacentMapping >::const_iterator it = mappings_.find( rank );
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
