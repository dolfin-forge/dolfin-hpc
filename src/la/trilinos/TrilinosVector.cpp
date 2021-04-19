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
  : Variable( "x", "a sparse vector" )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

Vector::Vector( size_t N, bool distributed )
  : Variable( "x", "a sparse vector" )
{
  // Create Tpetra vector
  init( N, distributed );
}

//-----------------------------------------------------------------------------

// Vector::Vector( TPVector x )
//   : Variable( "x", "a vector" )
// {
// }

//-----------------------------------------------------------------------------

Vector::Vector( Vector const & copy )
  : Variable( "x", "a vector" )
{
  if ( not copy.vec_.is_null() )
  {
    // Create with same map
    Teuchos::RCP< TPMap const > v_ghostmap( copy.vec_->getMap() );
    Teuchos::RCP< TPMap const > v_xmap( copy.vec_local_->getMap() );
    vec_ = Teuchos::rcp( new TPVector( v_ghostmap, 1 ) );

    vec_->assign( *copy.vec_ );
    vec_local_ = vec_->offsetViewNonConst( v_xmap, 0 );
  }
}

//-----------------------------------------------------------------------------

Vector::~Vector()
{
  clear();
}

//-----------------------------------------------------------------------------

auto Vector::copy() const -> Vector *
{
  return new Vector( *this );
}

//-----------------------------------------------------------------------------

auto Vector::zero() -> void
{
  dolfin_assert( not vec_.is_null() );
  vec_->putScalar( 0.0 );
}

//-----------------------------------------------------------------------------

auto Vector::apply( FinalizeType ) -> void
{
  dolfin_assert( not vec_.is_null() );

  // update_ghost_values();

  Teuchos::RCP< TPMap const > map = vec_local_->getMap();
  Teuchos::RCP< TPVector >    y( new TPVector( map, 1 ) );
  Teuchos::RCP< TPMap const > ghostmap = vec_->getMap();

  // Export from overlapping map ghostmap, to non-overlapping xmap
  Tpetra::Export< LO, GO, TPNode > exporter( ghostmap, map );

  // Forward export to reduction vector
  y->doExport( *vec_, exporter, Tpetra::ADD );

  // Copy back into _x_ghosted
  Tpetra::Import< LO, GO, TPNode > importer( map, ghostmap );
  vec_->doImport( *y, importer, Tpetra::INSERT );
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
  Teuchos::RCP< const TPMap > xmap = vec_->getMap();

  std::stringstream ss;

  size_t const rank = MPI::rank();
  size_t const m    = MPI::size();

  if ( rank == 0 )
  {
    ss << xmap->description() << "\n"
       << "trilinos::vector"
       << "\n---";
    for ( size_t j = 0; j != m; ++j )
      ss << "-";
    ss << "\n";
  }

  ss << rank << "] ";
  for ( size_t j = 0; j < m; ++j )
  {
    if ( xmap->isNodeGlobalElement( j ) )
      ss << "X";
    else
      ss << " ";
  }
  ss << "\n";

  Teuchos::ArrayRCP< real const > arr = vec_->getData( 0 );

  for ( std::size_t j = 0; j != xmap->getNodeNumElements(); ++j )
    ss << j << " -> " << xmap->getGlobalElement( j ) << " = " << arr[j] << "\n";
  ss << "\n";

  MPI_Comm mpi_comm = *( comm_->getRawMpiComm() );

  std::vector< std::string > out_str;
  MPIgather( mpi_comm, ss.str(), out_str );

  if ( rank == 0 )
  {
    for ( auto & s : out_str )
      std::cout << s;
  }

  // Teuchos::RCP< std::ostream > ostr;
  // Teuchos::FancyOStream os( ostr );
  // vec_->describe( os, Teuchos::VERB_EXTREME );

  vec_->print( std::cout );
}

//-----------------------------------------------------------------------------

auto Vector::size() const -> size_t
{
  if ( vec_.is_null() )
    return 0;
  else
    return vec_->getMap()->getMaxAllGlobalIndex() + 1;
}

//-----------------------------------------------------------------------------

auto Vector::local_size() const -> size_t
{
  if ( vec_.is_null() )
    return 0;
  else
    return vec_->getLocalLength();
}

//-----------------------------------------------------------------------------

auto Vector::offset() const -> size_t
{
  if ( vec_.is_null() )
    return 0;
  else
    return vec_->getMap()->getMinGlobalIndex();
}

//-----------------------------------------------------------------------------

auto Vector::init( GenericSparsityPattern const & sparsity_pattern ) -> void
{
  if ( vec_.is_null() )
  {
    SparsityPattern const & spattern = reinterpret_cast< SparsityPattern const & >( sparsity_pattern );

    size_t const pe_rank     = MPI::rank();
    size_t const nLocalRows  = spattern.size( 0 );
    size_t       nGlobalRows = 0;

    // compute the global vector size
    MPI::all_reduce< MPI::sum >( nLocalRows, nGlobalRows );
    dolfin_assert( nLocalRows <= nGlobalRows );

    // the local vector range is just our local range from spattern
    std::vector< GO > lIndices( nLocalRows );
    {
      std::vector< size_t > lRowRange( 2, DOLFIN_SIZE_T_MAX );
      spattern.get_range( pe_rank, lRowRange.data() );
      std::iota( lIndices.begin(), lIndices.end(), lRowRange[0] );
    }

    // global indices for our vector range map
    std::vector< GO > gIndices;

    if ( spattern.rank() == 1 )
    {
      // for rank == 1 the global indices are the same as our local indices
      gIndices = lIndices;
      // _ordered_set< size_t > gColset;
      // gColset.insert( lIndices.begin(), lIndices.end() );

      // std::vector< size_t > diag;
      // spattern.diagonal_entries( 0, diag );
      // gColset.insert( diag.begin(), diag.end() );

      // gIndices.resize( gColset.size() );
      // std::copy( gColset.begin(), gColset.end(), gIndices.begin() );
    }
    else // ( spattern.rank() == 2 )
    {
      _ordered_set< size_t > gColset;
      gColset.insert( lIndices.begin(), lIndices.end() );

      // dont be confused here, we take all the (unique) column indices from the
      // matrix to make sure we get all the ghost points into our vector range
      for ( size_t i = 0; i < nLocalRows; ++i )
      {
        // get columns for this row
        std::vector< size_t > diag, off_diag;
        spattern.diagonal_entries( i, diag );
        spattern.off_diagonal_entries( i, off_diag );

        // insert (potentially new) columns in the set
        gColset.insert( diag.begin(), diag.end() );
        gColset.insert( off_diag.begin(), off_diag.end() );
      }

      // copy the unique set of column indices into the column vector
      gIndices.resize( gColset.size() );
      std::copy( gColset.begin(), gColset.end(), gIndices.begin() );
    }

    Teuchos::ArrayView< GO > gColIndices_view( gIndices.data(), gIndices.size() );
    Teuchos::RCP< TPMap > map_global( new TPMap( nGlobalRows, gColIndices_view, indexBase, comm_ ) );

    // Vector - create with overlap
    vec_ = Teuchos::rcp( new TPVector( map_global, 1, true ) );

    // make sure we actually got a non-empty vector
    dolfin_assert( not vec_.is_null() );

    Teuchos::ArrayView< GO > lIndices_view( lIndices.data(), lIndices.size() );
    Teuchos::RCP< TPMap > map_local( new TPMap( nGlobalRows, lIndices_view, indexBase, comm_ ) );

    // Get a modifiable view into the ghosted vector
    vec_local_ = vec_->offsetViewNonConst( map_local, 0 );
  }
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N ) -> void
{
  init( N, true );
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N, bool ) -> void
{
  if ( vec_.is_null() )
  {
    // // compute the global vector size
    // size_t Nglobal = 0;
    // MPI::all_reduce< MPI::sum >( N, Nglobal );
    // dolfin_assert( N <= Nglobal );

    // Teuchos::RCP< TPMap > map( new TPMap( Nglobal, N, indexBase, comm_ ) );

    // // Vector - create with overlap
    // vec_= Teuchos::rcp( new TPVector( map, 1, true ) );

    // // make sure we actually got a non-empty vector
    // dolfin_assert( not vec_.is_null() );

    // // Get a modifiable view into the ghosted vector
    // vec_local_ = vec_->offsetViewNonConst( map, 0 );

    size_t Nglobal = 0;
    MPI::all_reduce< MPI::sum >( N, Nglobal );
    dolfin_assert( N <= Nglobal );

    std::vector< size_t > Nall( MPI::size() );
    MPI::all_gather( &N, 1, Nall.data(), 1 );

    size_t my_start = 0;
    for ( size_t i = 0; i < Nall.size(); ++i )
      if ( i < MPI::rank() )
        my_start += Nall[i];

    std::vector< size_t > global( Nglobal );
    std::iota( global.begin(), global.end(), 0 );

    std::vector< size_t > local( N );
    std::iota( local.begin(), local.end(), my_start );

    Teuchos::ArrayView< GO > g_view( global.data(), global.size() );
    Teuchos::RCP< TPMap > map_global( new TPMap( Nglobal, g_view, indexBase, comm_ ) );

    // Vector - create with overlap
    vec_ = Teuchos::rcp( new TPVector( map_global, 1, true ) );

    // make sure we actually got a non-empty vector
    dolfin_assert( not vec_.is_null() );

    Teuchos::ArrayView< GO > l_view( local.data(), local.size() );
    Teuchos::RCP< TPMap > map_local( new TPMap( Nglobal, l_view, indexBase, comm_ ) );

    // Get a modifiable view into the ghosted vector
    vec_local_ = vec_->offsetViewNonConst( map_local, 0 );
  }
}

//-----------------------------------------------------------------------------

auto Vector::init_ghosted( size_t                           n,
                           _ordered_set< size_t > &         indices,
                           _ordered_map< size_t, size_t > & map ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Vector::get( real * values ) const -> void
{
  dolfin_assert( not vec_.is_null() );

  message( "get( values )" );
  dolfin_assert( values != nullptr );

  Teuchos::ArrayRCP< const real > arr = vec_->getData( 0 );
  std::copy( arr.get(), arr.get() + local_size(), values );
}

//-----------------------------------------------------------------------------

auto Vector::set( real * values ) -> void
{
  size_t const num_values = local_size();

  message( "set( values )" );

  if ( num_values != 0 )
  {
    dolfin_assert( not vec_.is_null() );
    dolfin_assert( values != nullptr );

    Teuchos::ArrayRCP< real > arr = vec_->getDataNonConst( 0 );
    std::copy( values, values + num_values, arr.get() );
  }
}

//-----------------------------------------------------------------------------

auto Vector::add( real * values ) -> void
{
  dolfin_assert( not vec_.is_null() );

  message( "add( values )" );

  size_t const num_values = local_size();

  for ( size_t i = 0; i < num_values; ++i )
    vec_->sumIntoLocalValue( i, 0, values[i] );
}

//-----------------------------------------------------------------------------

auto Vector::get( real * block, size_t m, const size_t * rows ) const -> void
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::RCP< TPMap const >     map = vec_->getMap();
  Teuchos::ArrayRCP< real const > arr = vec_->getData( 0 );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( not map->isNodeGlobalElement( rows[i] ) )
    {
      std::string msg = "trilinos::Vector::get(): Row " + std::to_string(rows[i]);
      msg += " is not global (process " + std::to_string( MPI::rank() ) + ")";
      // warning( "trilinos::Vector::get(): Row %d is not valid", rows[i] );
      std::cout << msg << std::endl;
    }

    LO const idx = map->getLocalElement( rows[i] );
    if ( idx != Teuchos::OrdinalTraits< LO >::invalid() )
    {
      block[i] = arr[idx];
    }
    else
    {
      std::string msg = "trilinos::Vector::get(): Row " + std::to_string(rows[i]);
      msg += " is not local (process " + std::to_string( MPI::rank() ) + ")";
      // warning( "trilinos::Vector::get(): Row %d is not valid", rows[i] );
      std::cout << msg << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::set( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::RCP< TPMap const > map = vec_->getMap();

  for ( size_t i = 0; i < m; ++i )
  {
    if ( map->isNodeGlobalElement( rows[i] ) )
    {
      vec_->replaceGlobalValue( rows[i], 0, block[i] );
    }
    else
    {
      std::string msg = "trilinos::Vector::set(): Row " + std::to_string(rows[i]);
      msg += " is not local (process " + std::to_string( MPI::rank() ) + ")";
      // warning( "trilinos::Vector::set(): Row %d is not valid", rows[i] );
      std::cout << msg << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::add( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::RCP< TPMap const > map = vec_->getMap();

  for ( size_t i = 0; i < m; ++i )
  {
    if ( map->isNodeGlobalElement( rows[i] ) )
    {
      vec_->sumIntoGlobalValue( rows[i], 0, block[i] );
    }
    else
    {
      std::string msg = "trilinos::Vector::add(): Row " + std::to_string(rows[i]);
      msg += " is not local (process " + std::to_string( MPI::rank() ) + ")";
      // warning( "trilinos::Vector::add(): Row %d is not local", rows[i] );
      std::cout << msg << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::axpy( real a, const GenericVector & x ) -> void
{
  this->axpby( a, x, 1.0 );
}

//-----------------------------------------------------------------------------

auto Vector::axpby( real a, const GenericVector & x, real b ) -> void
{
  dolfin_assert( not vec_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.vec_.is_null() );

  vec_->update( a, *X_.vec_, b );
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
  dolfin_assert( not vec_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.vec_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_.is_null() );

  vec_->update( a, *X_.vec_, b, *Y_.vec_, c );
}

//-----------------------------------------------------------------------------

auto Vector::inner( const GenericVector & y ) const -> real
{
  dolfin_assert( not vec_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_.is_null() );

  std::vector< real >              val( 1 );
  const Teuchos::ArrayView< real > result( val );

  vec_->dot( *Y_.vec_, result );

  return val[0];
}

//-----------------------------------------------------------------------------

auto Vector::norm( VectorNormType type ) const -> real
{
  dolfin_assert( not vec_.is_null() );
  using TPMagType = Tpetra::Vector<>::mag_type;

  std::vector< TPMagType >              norms( 1 );
  Teuchos::ArrayView< TPMagType > const norm_view( norms );

  switch ( type )
  {
    case l1:
      vec_->norm1( norm_view );
      break;
    case l2:
      vec_->norm2( norm_view );
      break;
    default:
      vec_->normInf( norm_view );
      break;
  }

  return norms[0];
}

//-----------------------------------------------------------------------------

auto Vector::min() const -> real
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::ArrayRCP< real const > arr = vec_->getData( 0 );
  real min = *std::min_element( arr.get(), arr.get() + arr.size() );

  MPI::all_reduce_in_place< MPI::min >( min );

  return min;
}

//-----------------------------------------------------------------------------

auto Vector::max() const -> real
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::ArrayRCP< real const > arr = vec_->getData( 0 );
  real max = *std::min_element( arr.get(), arr.get() + arr.size() );

  MPI::all_reduce_in_place< MPI::max >( max );

  return max;
}

//-----------------------------------------------------------------------------

auto Vector::pointwise( const GenericVector & x,
                        VectorPointwiseOp     op ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Vector::operator*=( const real a ) -> Vector &
{
  dolfin_assert( not vec_.is_null() );

  vec_->scale( a );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator/=( const real a ) -> Vector &
{
  dolfin_assert( not vec_.is_null() );
  dolfin_assert( a != 0.0 );

  ( *this ) *= 1.0 / a;

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator*=( const GenericVector & y ) -> Vector &
{
  dolfin_assert( not vec_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.vec_.is_null() );

  vec_->elementWiseMultiply( 1.0, *( vec_->getVector( 0 ) ), *( Y_.vec_ ), 0.0 );

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
    dolfin_assert( not v.vec_.is_null() );
    dolfin_assert( not vec_.is_null() );

    vec_->assign( *v.vec_ );
  }

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator=( real a ) -> Vector &
{
  dolfin_assert( not vec_.is_null() );

  vec_->putScalar( a );

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

void Vector::update_ghost_values()
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::RCP< TPMap const > localmap( vec_local_->getMap() );
  Teuchos::RCP< TPMap const > map( vec_->getMap() );

  // Export from non-overlapping map x, to overlapping ghost values
  Tpetra::Import< LO, GO, TPNode > importer( localmap, map );

  // FIXME: is this safe, since vec_local_ is a view into vec_?
  vec_->doImport( *vec_local_, importer, Tpetra::INSERT );
}

//-----------------------------------------------------------------------------

auto Vector::clear() -> void
{
  vec_ = Teuchos::null;
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
