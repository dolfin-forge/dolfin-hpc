// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosVector.h>

#include <dolfin/la/trilinos/TrilinosFactory.h>

#include <numeric>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

Vector::Vector()
  : Variable( "x", "a vector" )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

Vector::Vector( size_t N, bool distributed )
  : Variable( "x", "a vector" )
{
  // Create Tpetra vector
  init( N, distributed );
}

//-----------------------------------------------------------------------------

Vector::Vector( Vector const & copy )
  : Variable( "x", "a vector" )
{
  if ( not copy.vec_local_.is_null() )
  {
    // Create with same map
    Teuchos::RCP< TPMap const > v_local( copy.vec_local_->getMap() );
    vec_local_ = Teuchos::rcp( new TPVector( v_local, 1 ) );
    vec_local_->assign( *copy.vec_local_ );
  }

  if ( not copy.vec_ghost_.is_null() )
  {
    Teuchos::RCP< TPMap const > v_ghost( copy.vec_ghost_->getMap() );
    vec_ghost_ = Teuchos::rcp( new TPVector( v_ghost, 1 ) );
    vec_ghost_->assign( *copy.vec_ghost_ );
  }
}

//-----------------------------------------------------------------------------

Vector::~Vector()
{
  vec_local_ = Teuchos::null;
  vec_ghost_ = Teuchos::null;
}

//-----------------------------------------------------------------------------

auto Vector::copy() const -> Vector *
{
  return new Vector( *this );
}

//-----------------------------------------------------------------------------

auto Vector::zero() -> void
{
  *this = 0.0;
}

//-----------------------------------------------------------------------------

auto Vector::apply( FinalizeType ) -> void
{
  if ( not vec_ghost_.is_null() )
  {
    dolfin_assert( not vec_local_.is_null() );

    Teuchos::RCP< TPMap const > map_ghost( vec_ghost_->getMap() );
    Teuchos::RCP< TPMap const > map_local( vec_local_->getMap() );

    // check if ghosts have been modified on any process
    bool ghosts_need_sync = not ghost_stash_.empty();
    MPI::all_reduce_in_place< MPI::sum >( ghosts_need_sync, MPI::DOLFIN_COMM );

    if ( ghosts_need_sync )
    {
      // set ghosts to zero
      vec_ghost_->putScalar( 0.0 );

      // insert ghost additions
      for ( std::pair< GO const, real > const & ghost : ghost_stash_ )
        vec_ghost_->replaceGlobalValue( ghost.first, 0, ghost.second );

      // clear ghost stash
      ghost_stash_.clear();

      // add our local ghost part to the global vector
      Tpetra::Export< LO, GO, TPNode > exporter( map_ghost, map_local );
      vec_local_->doExport( *vec_ghost_, exporter, Tpetra::ADD );
    }

    // replace our ghost points by the global owners value
    Tpetra::Export< LO, GO, TPNode > exporter( map_local, map_ghost );
    vec_ghost_->doExport( *vec_local_, exporter, Tpetra::REPLACE );
  }
}

//-----------------------------------------------------------------------------

auto MPIgather( MPI_Comm                     comm,
                const std::string &          in_values,
                std::vector< std::string > & out_values,
                unsigned int                 receiving_process = 0 ) -> void
{
#ifdef DOLFIN_HAVE_MPI
  size_t const comm_size = MPI::size();

  // Get data size on each process
  std::vector< int > pcounts( comm_size );
  int local_size = in_values.size();

  MPI_Gather( &local_size, 1, MPI_type< int >::value,
              pcounts.data(), 1, MPI_type< int >::value,
              receiving_process, comm );

  // Build offsets
  std::vector< int > offsets( comm_size + 1, 0 );
  for ( std::size_t i = 1; i <= comm_size; ++i )
    offsets[i] = offsets[i - 1] + pcounts[i - 1];

  // Gather
  size_t const n = std::accumulate( pcounts.begin(), pcounts.end(), 0 );
  std::vector< char > _out( n );

  MPI_Gatherv( const_cast< char * >( in_values.data() ),
               in_values.size(), MPI_type< char >::value, _out.data(),
               pcounts.data(), offsets.data(), MPI_type< char >::value,
               receiving_process, comm );

  // Rebuild
  out_values.resize( comm_size );
  for ( std::size_t p = 0; p < comm_size; ++p )
  {
    out_values[p] = std::string( _out.begin() + offsets[p],
                                 _out.begin() + offsets[p + 1] );
  }
#else
  out_values.clear();
  out_values.push_back( in_values );
#endif
}

auto Vector::disp( size_t ) const -> void
{
  std::stringstream ss;

  size_t const rank = MPI::rank();
  size_t const m    = MPI::size();

  if ( not vec_local_.is_null() )
  {
    Teuchos::RCP< const TPMap > map_local = vec_local_->getMap();

    if ( rank == 0 )
    {
      ss << map_local->description() << "\n"
         << "trilinos::vector"
         << "\n---";
      for ( size_t j = 0; j < m; ++j )
        ss << "-";
      ss << "\n";
    }

    ss << rank << "] local:\n";
    ss << map_local->description() << "\n";
    Teuchos::ArrayRCP< real const > arr_local = vec_local_->getData( 0 );

    for ( std::size_t j = 0; j != map_local->getNodeNumElements(); ++j )
      ss << j << " -> " << map_local->getGlobalElement( j ) << " = " << arr_local[j] << "\n";
    ss << "\n";
  }

  if ( not vec_ghost_.is_null() )
  {
    ss << rank << "] ghost:\n";
    Teuchos::RCP< const TPMap > map_ghost     = vec_ghost_->getMap();
    Teuchos::ArrayRCP< real const > arr_ghost = vec_ghost_->getData( 0 );
    ss << map_ghost->description() << "\n";

    for ( std::size_t j = 0; j != map_ghost->getNodeNumElements(); ++j )
      ss << j << " -> " << map_ghost->getGlobalElement( j ) << " = " << arr_ghost[j] << "\n";
    ss << "\n";
  }

  MPI_Comm mpi_comm = *( comm_->getRawMpiComm() );

  std::vector< std::string > out_str;
  MPIgather( mpi_comm, ss.str(), out_str );

  if ( rank == 0 )
  {
    for ( auto & s : out_str )
      std::cout << s;
  }

  // FIXME use just this?
  // Teuchos::RCP< std::ostream > ostr;
  // Teuchos::FancyOStream os( ostr );
  // vec_->describe( os, Teuchos::VERB_EXTREME );

  vec_local_->print( std::cout );
}

//-----------------------------------------------------------------------------

auto Vector::size() const -> size_t
{
  if ( vec_local_.is_null() )
    return 0;
  else
    return vec_local_->getMap()->getMaxAllGlobalIndex() + 1;
}

//-----------------------------------------------------------------------------

auto Vector::local_size() const -> size_t
{
  if ( vec_local_.is_null() )
    return 0;
  else
    return vec_local_->getLocalLength();
}

//-----------------------------------------------------------------------------

auto Vector::offset() const -> size_t
{
  if ( vec_local_.is_null() )
    return 0;
  else
    return vec_local_->getMap()->getMinGlobalIndex();
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N ) -> void
{
  init( N, true );
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N, bool ) -> void
{
  if ( vec_local_.is_null() )
  {
    // compute global number of elements amd local offset
    std::vector< size_t > Nall( MPI::size() );
    MPI::all_gather( &N, 1, Nall.data(), 1 );

    size_t const Nglobal = std::accumulate( Nall.begin(), Nall.end(), 0 );
    size_t const offset  = std::accumulate( Nall.data(), Nall.data() + MPI::rank(), 0 );

    dolfin_assert( N <= Nglobal );
    dolfin_assert( offset <= Nglobal );

    // create our local index range
    std::vector< size_t > indices( N );
    std::iota( indices.begin(), indices.end(), offset );

    // local and global map are the same in this case
    Teuchos::ArrayView< GO > indices_view( indices.data(), indices.size() );
    Teuchos::RCP< TPMap > map( new TPMap( Nglobal, indices_view, indexBase, comm_ ) );

    // Vector - create with overlap
    vec_local_ = Teuchos::rcp( new TPVector( map, 1, true ) );

    // make sure we actually got a non-empty vector
    dolfin_assert( not vec_local_.is_null() );

    // for now we don't have any ghost info
    vec_ghost_ = Teuchos::null;
  }
}

//-----------------------------------------------------------------------------

auto Vector::init_ghosted( _ordered_set< size_t > & indices ) -> void
{
  dolfin_assert( not vec_local_.is_null() );

  // sync ghosts -> local
  if ( not vec_ghost_.is_null() )
    apply();

  // compute the ghost map ( global map - local map )
  std::vector< GO > new_ghost_indices;
  {
    for ( size_t idx : indices )
    {
      if ( not vec_local_->getMap()->isNodeGlobalElement( idx ) )
      {
        new_ghost_indices.push_back( idx );
      }
    }
  }
  Teuchos::ArrayView< GO > new_ghost_indices_view( new_ghost_indices.data(),
                                                   new_ghost_indices.size() );
  Teuchos::RCP< TPMap > map_ghost( new TPMap( this->size(), new_ghost_indices_view,
                                              indexBase, comm_ ) );

  // set the ghost vector
  vec_ghost_ = Teuchos::rcp( new TPVector( map_ghost, 1, true ) );

  // sync local -> ghost
  apply();
}

//-----------------------------------------------------------------------------

auto Vector::get( real * block, size_t m, const size_t * rows ) const -> void
{
  dolfin_assert( not vec_local_.is_null() );

  Teuchos::ArrayRCP< real const > local_data = vec_local_->getData( 0 );
  Teuchos::ArrayRCP< real const > ghost_data = ( not vec_ghost_.is_null() )
                                               ? vec_ghost_->getData( 0 )
                                               : Teuchos::null;

  for ( size_t i = 0; i < m; ++i )
  {
    if ( vec_local_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      LO const idx = vec_local_->getMap()->getLocalElement( rows[i] );
      block[i]     = local_data[idx];
    }
    else if ( not vec_ghost_.is_null() )
    {
      if ( vec_ghost_->getMap()->isNodeGlobalElement( rows[i] ) )
      {
        LO const idx = vec_ghost_->getMap()->getLocalElement( rows[i] );
        block[i]     = ghost_data[idx];
      }
    }
    else
    {
      error( "trilinos::Vector::get(): invalid row id: %d", rows[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::set( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_local_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( vec_local_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      vec_local_->replaceGlobalValue( rows[i], 0, block[i] );
    }
    else if ( not vec_ghost_.is_null() )
    {
      if( vec_ghost_->getMap()->isNodeGlobalElement( rows[i] ) )
      {
        ghost_stash_[rows[i]] = block[i];
      }
    }
    else
    {
      error( "trilinos::Vector::set(): invalid row id: %d (value: %f)",
             rows[i], block[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::add( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_local_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( vec_local_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      vec_local_->sumIntoGlobalValue( rows[i], 0, block[i] );
    }
    else if ( not vec_ghost_.is_null() )
    {
      if( vec_ghost_->getMap()->isNodeGlobalElement( rows[i] ) )
      {
        ghost_stash_[rows[i]] += block[i];
      }
    }
    else
    {
      error( "trilinos::Vector::add(): invalid row id: %d (value: %f)",
             rows[i], block[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::get( real * values ) const -> void
{
  dolfin_assert( not vec_local_.is_null() );
  dolfin_assert( values != nullptr );

  Teuchos::ArrayRCP< real const > data = vec_local_->getData( 0 );
  std::copy( data.get(), data.get() + local_size(), values );
}

//-----------------------------------------------------------------------------

auto Vector::set( real * values ) -> void
{
  dolfin_assert( not vec_local_.is_null() );
  dolfin_assert( values != nullptr );

  size_t const num_values = local_size();

  for ( size_t i = 0; i < num_values; ++i )
    vec_local_->replaceLocalValue( i, 0, values[i] );
}

//-----------------------------------------------------------------------------

auto Vector::add( real * values ) -> void
{
  dolfin_assert( not vec_local_.is_null() );
  dolfin_assert( values != nullptr );

  size_t const num_values = local_size();

  for ( size_t i = 0; i < num_values; ++i )
    vec_local_->sumIntoLocalValue( i, 0, values[i] );
}

//-----------------------------------------------------------------------------

auto Vector::axpy( real a, const GenericVector & x ) -> void
{
  this->axpby( a, x, 1.0 );
}

//-----------------------------------------------------------------------------

auto Vector::axpby( real a, const GenericVector & x, real b ) -> void
{
  dolfin_assert( not vec_local_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.vec_local_.is_null() );

  vec_local_->update( a, *X_.vec_local_, b );
}

//-----------------------------------------------------------------------------

auto Vector::waxpy( real a, const GenericVector & x, const GenericVector & y ) -> void
{
  this->axpbypcz( a, x, 1.0, y, 0.0 );
}

//-----------------------------------------------------------------------------

auto Vector::axpbypcz( real                  a,
                       const GenericVector & x,
                       real                  b,
                       const GenericVector & y,
                       real                  c ) -> void
{
  dolfin_assert( not vec_local_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.vec_local_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_local_.is_null() );

  vec_local_->update( a, *X_.vec_local_, b, *Y_.vec_local_, c );
}

//-----------------------------------------------------------------------------

auto Vector::inner( const GenericVector & y ) const -> real
{
  dolfin_assert( not vec_local_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_local_.is_null() );

  std::vector< real >              value( 1 );
  Teuchos::ArrayView< real > const result( value );

  vec_local_->dot( *Y_.vec_local_, result );

  return value[0];
}

//-----------------------------------------------------------------------------

auto Vector::norm( VectorNormType type ) const -> real
{
  dolfin_assert( not vec_local_.is_null() );
  using TPMagType = Tpetra::Vector<>::mag_type;

  std::vector< TPMagType >              norm( 1 );
  Teuchos::ArrayView< TPMagType > const norm_view( norm );

  switch ( type )
  {
    case l1:
      vec_local_->norm1( norm_view );
      break;
    case l2:
      vec_local_->norm2( norm_view );
      break;
    default:
      vec_local_->normInf( norm_view );
      break;
  }

  return norm[0];
}

//-----------------------------------------------------------------------------

auto Vector::min() const -> real
{
  dolfin_assert( not vec_local_.is_null() );

  Teuchos::ArrayRCP< real const > data = vec_local_->getData( 0 );
  real min = *std::min_element( data.get(), data.get() + data.size() );

  MPI::all_reduce_in_place< MPI::min >( min );

  return min;
}

//-----------------------------------------------------------------------------

auto Vector::max() const -> real
{
  dolfin_assert( not vec_local_.is_null() );

  Teuchos::ArrayRCP< real const > data = vec_local_->getData( 0 );
  real max = *std::max_element( data.get(), data.get() + data.size() );

  MPI::all_reduce_in_place< MPI::max >( max );

  return max;
}

//-----------------------------------------------------------------------------

auto Vector::pointwise( const GenericVector &, VectorPointwiseOp ) const -> void
{
  // FIXME
  // elementWiseMultiply  ?!
}

//-----------------------------------------------------------------------------

auto Vector::operator*=( const real a ) -> Vector &
{
  if ( not vec_local_.is_null() )
    vec_local_->scale( a );

  if ( not vec_ghost_.is_null() )
    vec_ghost_->scale( a );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator/=( const real a ) -> Vector &
{
  dolfin_assert( not vec_local_.is_null() );
  dolfin_assert( a != 0.0 );

  ( *this ) *= 1.0 / a;

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator*=( const GenericVector & y ) -> Vector &
{
  dolfin_assert( not vec_local_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_local_.is_null() );

  if ( not vec_local_.is_null() )
    vec_local_->elementWiseMultiply( 1.0, *( vec_local_->getVector( 0 ) ),
                                     *( Y_.vec_local_ ), 0.0 );

  if ( not vec_ghost_.is_null() and not Y_.vec_ghost_.is_null() )
    vec_ghost_->elementWiseMultiply( 1.0, *( vec_ghost_->getVector( 0 ) ),
                                     *( Y_.vec_ghost_ ), 0.0 );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator+=( const GenericVector & x ) -> Vector &
{
  this->axpy( 1.0, x );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator-=( const GenericVector & x ) -> Vector &
{
  this->axpy( -1.0, x );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator=( const GenericVector & v ) -> Vector &
{
  *this = v.down_cast< trilinos::Vector >();
  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator=( Vector const & v ) -> Vector &
{
  // Check that vector lengths are equal
  if ( size() != v.size() )
  {
    error( "trilinos::Vector: Vectors must be of equal size for assignment" );
  }

  // Check that vector local ranges are equal (relevant in parallel)
  if ( local_size() != v.local_size() )
  {
    error( "trilinos::Vector: Vectors must equal local size for assignment" );
  }

  // Check for self-assignment
  if ( this != &v )
  {
    // Copy data (local operation)
    dolfin_assert( not v.vec_local_.is_null() );
    dolfin_assert( not vec_local_.is_null() );

    vec_local_->assign( *v.vec_local_ );

    if ( not vec_local_.is_null() )
      vec_local_->assign( *v.vec_local_ );

    if ( not vec_ghost_.is_null() and not v.vec_ghost_.is_null() )
      vec_ghost_->assign( *v.vec_ghost_ );
  }

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator=( real a ) -> Vector &
{
  if ( not vec_local_.is_null() )
    vec_local_->putScalar( a );

  if ( not vec_ghost_.is_null() )
    vec_ghost_->putScalar( a );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::vec() const -> TPVectorPtr
{
  return vec_local_;
}

//-----------------------------------------------------------------------------

auto Vector::factory() const -> LinearAlgebraFactory &
{
  return trilinos::Factory::instance();
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
