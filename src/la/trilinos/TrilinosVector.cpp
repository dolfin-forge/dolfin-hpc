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
  vec_ = Teuchos::rcp( new TPVector( copy.vec_->getMap(),
                                     copy.vec_->getNumVectors(),
                                     false ) );
  vec_->assign( *copy.vec_ );
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

  // FIXME
  // update_ghost_values();

  Teuchos::RCP< TPMap const > map = vec_->getMap();
  Teuchos::RCP< TPVector >    y( new TPVector( map, 1 ) );
  Teuchos::RCP< TPMap const > ghostmap = vec_->getMap();

  // Export from overlapping map ghostmap, to non-overlapping xmap
  Tpetra::Export< LO, GO, TPNode > exporter( ghostmap, map );

  // Forward export to reduction vector
  y->doExport( *vec_, exporter, Tpetra::ADD );

  // Copy back into _x_ghosted
  Tpetra::Import< LO, GO, TPNode > importer( map, map );
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

auto Vector::init( size_t N ) -> void
{
  init( N, true );
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N, bool ) -> void
{
  if ( vec_.is_null() )
  {
    // compute the global vector size
    size_t Nglobal = 0;
    MPI::all_reduce< MPI::sum >( N, Nglobal );
    dolfin_assert( N <= Nglobal );

    Teuchos::RCP< TPMap > map( new TPMap( Nglobal, N, indexBase, comm_ ) );

    // Vector - create with overlap
    vec_= Teuchos::rcp( new TPVector( map, 1, true ) );

    // make sure we actually got a non-empty vector
    dolfin_assert( not vec_.is_null() );
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
  dolfin_assert( values != nullptr );

  Teuchos::ArrayRCP< const real > arr = vec_->getData( 0 );
  std::copy( arr.get(), arr.get() + local_size(), values );
}

//-----------------------------------------------------------------------------

auto Vector::set( real * values ) -> void
{
  size_t const num_values = local_size();

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

  size_t const num_values = local_size();

  for ( size_t i = 0; i < num_values; ++i )
    vec_->sumIntoLocalValue( i, 0, values[i] );
}

//-----------------------------------------------------------------------------

auto Vector::get( real * block, size_t m, const size_t * rows ) const -> void
{
  dolfin_assert( not vec_.is_null() );

  Teuchos::RCP< TPMap const >     xmap = vec_->getMap();
  Teuchos::ArrayRCP< real  const> xarr = vec_->getData( 0 );

  for ( size_t i = 0; i < m; ++i )
  {
    LO const idx = xmap->getLocalElement( rows[i] );
    if ( idx != Teuchos::OrdinalTraits< LO >::invalid() )
    {
      block[i] = xarr[idx];
    }
    else
    {
      error( "trilinos::Vector: Row %d is not valid", rows[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::set( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( vec_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      vec_->replaceGlobalValue( rows[i], 0, block[i] );
    }
    else
    {
      error( "trilinos::Vector: Row %d is not valid", rows[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::add( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not vec_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( vec_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      vec_->sumIntoGlobalValue( rows[i], 0, block[i] );
    }
    else
    {
      error( "trilinos::Vector: Row %d is not local", rows[i] );
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
  return vec_;
}

//-----------------------------------------------------------------------------

auto Vector::factory() const -> LinearAlgebraFactory &
{
  return trilinos::Factory::instance();
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
