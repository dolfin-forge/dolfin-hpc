// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosVector.h>

#include <dolfin/la/trilinos/TrilinosFactory.h>

// from the trilinos bug tracker:
// "you shouldn't ever need to include a "_def""
// well... that was a lie...
#include <Tpetra_DistObject.hpp>
#include <Tpetra_DistObject_def.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Export_def.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Import_def.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_ImportExportData_def.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Vector_def.hpp>

#include <Tpetra_Details_Transfer.hpp>
#include <Tpetra_Details_Transfer_def.hpp>

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
  if ( not copy.x_.is_null() )
  {
    Teuchos::RCP< TPMap const > copy_ghostmap( copy.x_ghosted_->getMap() );
    Teuchos::RCP< TPMap const > copy_xmap( copy.x_->getMap() );

    x_ghosted_ = Teuchos::rcp( new TPVector( copy_ghostmap, 1 ) );
    x_ghosted_->assign( *copy.x_ghosted_ );

    x_ = x_ghosted_->offsetViewNonConst( copy_xmap, 0 );
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
  dolfin_assert( not x_ghosted_.is_null() );
  x_ghosted_->putScalar( 0.0 );
}

//-----------------------------------------------------------------------------

auto Vector::apply( FinalizeType ) -> void
{
  dolfin_assert( not x_.is_null() );

  // update_ghost_values();

  Teuchos::RCP< const TPMap > xmap = x_->getMap();
  Teuchos::RCP< TPVector >    y( new TPVector( xmap, 1 ) );
  Teuchos::RCP< const TPMap > ghostmap = x_ghosted_->getMap();

  // Export from overlapping map ghostmap, to non-overlapping xmap
  Tpetra::Export< TPVector::local_ordinal_type,
                  TPVector::global_ordinal_type,
                  TPVector::node_type > exporter( ghostmap, xmap );

  // Forward export to reduction vector
  y->doExport( *x_ghosted_, exporter, Tpetra::ADD );

  // Copy back into _x_ghosted
  Tpetra::Import< TPVector::local_ordinal_type,
                  TPVector::global_ordinal_type,
                  TPVector::node_type > importer( xmap, ghostmap );
  x_ghosted_->doImport( *y, importer, Tpetra::INSERT );
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
  Teuchos::RCP< const TPMap > xmap = x_->getMap();

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

  for ( std::size_t j = 0; j != xmap->getNodeNumElements(); ++j )
    ss << j << " -> " << xmap->getGlobalElement( j ) << "\n";
  ss << "\n";

  const Teuchos::RCP< const Teuchos::MpiComm< int > > _mpi_comm =
    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm< int > >(
      xmap()->getComm() );
  MPI_Comm mpi_comm = *( _mpi_comm->getRawMpiComm() );

  std::vector< std::string > out_str;
  MPIgather( mpi_comm, ss.str(), out_str );

  if ( rank == 0 )
  {
    for ( auto & s : out_str )
      std::cout << s;
  }
}

//-----------------------------------------------------------------------------

auto Vector::size() const -> size_t
{
  if ( x_.is_null() )
    return 0;
  else
    return x_->getMap()->getMaxAllGlobalIndex() + 1;
}

//-----------------------------------------------------------------------------

auto Vector::local_size() const -> size_t
{
  if ( x_.is_null() )
    return 0;
  else
    return x_->getLocalLength();
}

//-----------------------------------------------------------------------------

auto Vector::offset() const -> size_t
{
  if ( x_.is_null() )
    return 0;
  else
    return x_->getMap()->getMinGlobalIndex();
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N ) -> void
{
  init( N, true );
}

//-----------------------------------------------------------------------------

auto Vector::init( size_t N, bool distributed ) -> void
{
  if ( not distributed )
  {
    // FIXME this should also work with just 1 process
    error( "trilinos::Vector: only implemented for distributed scenarios" );
  }

  if ( not x_.is_null() )
  {
    error( "trilinos::Vector: Vector cannot be initialised more than once" );
  }

  size_t Nglobal = 0;
  MPI::all_reduce< MPI::sum >( N, Nglobal );
  dolfin_assert( N <= Nglobal );

  // FIXME use the correct comm here
  Teuchos::RCP< const Teuchos::MpiComm< int > > _comm(
    new Teuchos::MpiComm< int >( Teuchos::MpiComm< int >( MPI::DOLFIN_COMM ) ) );

  Teuchos::RCP< TPMap > _map( new TPMap( Nglobal, N, indexBase, _comm ) );
  Teuchos::RCP< TPMap > _ghost_map;
  std::vector< GO > local_to_global_map;

  // Save a map for the ghosting of values on other processes
  if ( local_to_global_map.size() != 0 )
  {
    const Teuchos::ArrayView< GO const > local_indices( local_to_global_map );
    _ghost_map = Teuchos::rcp( new TPMap( N, local_indices, indexBase, _comm ) );
  }
  else
  {
    _ghost_map = _map;
  }

  // Vector - create with overlap
  x_ghosted_ = Teuchos::rcp( new TPVector( _ghost_map, 1 ) );

  // make sure we actually got a non-empty vector
  dolfin_assert( not x_ghosted_.is_null() );

  // Get a modifiable view into the ghosted vector
  x_ = x_ghosted_->offsetViewNonConst( _map, 0 );
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
  dolfin_assert( not x_.is_null() );
  dolfin_assert( values != nullptr );

  Teuchos::ArrayRCP< const real > arr = x_->getData( 0 );
  std::copy( arr.get(), arr.get() + local_size(), values );
}

//-----------------------------------------------------------------------------

auto Vector::set( real * values ) -> void
{
  size_t const num_values = local_size();

  if ( num_values != 0 )
  {
    dolfin_assert( not x_.is_null() );
    dolfin_assert( values != nullptr );

    Teuchos::ArrayRCP< real > arr = x_->getDataNonConst( 0 );
    std::copy( values, values + num_values, arr.get() );
  }
}

//-----------------------------------------------------------------------------

auto Vector::add( real * values ) -> void
{
  dolfin_assert( not x_.is_null() );

  size_t const num_values = local_size();

  for ( size_t i = 0; i < num_values; ++i )
    x_->sumIntoLocalValue( i, 0, values[i] );
}

//-----------------------------------------------------------------------------

auto Vector::get( real * block, size_t m, const size_t * rows ) const -> void
{
  dolfin_assert( not x_ghosted_.is_null() );

  Teuchos::RCP< const TPMap >     xmap = x_ghosted_->getMap();
  Teuchos::ArrayRCP< const real > xarr = x_ghosted_->getData( 0 );

  for ( size_t i = 0; i < m; ++i )
  {
    const int idx = xmap->getLocalElement( rows[i] );
    if ( idx != Teuchos::OrdinalTraits< int >::invalid() )
      block[i] = xarr[idx];
    else
    {
      error( "trilinos::Vector: Row %d is not valid", rows[i] );
    }
  }
}

//-----------------------------------------------------------------------------

auto Vector::set( const real * block, size_t m, const size_t * rows ) -> void
{
  dolfin_assert( not x_ghosted_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( x_ghosted_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      x_ghosted_->replaceGlobalValue( rows[i], 0, block[i] );
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
  dolfin_assert( not x_ghosted_.is_null() );

  for ( size_t i = 0; i < m; ++i )
  {
    if ( x_ghosted_->getMap()->isNodeGlobalElement( rows[i] ) )
    {
      x_ghosted_->sumIntoGlobalValue( rows[i], 0, block[i] );
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
  dolfin_assert( not x_ghosted_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.x_ghosted_.is_null() );

  x_ghosted_->update( a, *X_.x_ghosted_, b );
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
  dolfin_assert( not x_ghosted_.is_null() );

  trilinos::Vector const & X_ = x.down_cast< trilinos::Vector >();
  dolfin_assert( not X_.x_ghosted_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.x_ghosted_.is_null() );

  x_ghosted_->update( a, *X_.x_ghosted_, b, *Y_.x_ghosted_, c );
}

//-----------------------------------------------------------------------------

auto Vector::inner( const GenericVector & y ) const -> real
{
  dolfin_assert( not x_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.x_ghosted_.is_null() );

  std::vector< real >              val( 1 );
  const Teuchos::ArrayView< real > result( val );

  x_->dot( *Y_.x_, result );

  return val[0];
}

//-----------------------------------------------------------------------------

auto Vector::norm( VectorNormType type ) const -> real
{
  dolfin_assert( not x_.is_null() );
  using TPMagType = Tpetra::MultiVector<>::mag_type;

  std::vector< TPMagType >              norms( 1 );
  Teuchos::ArrayView< TPMagType > const norm_view( norms );

  switch ( type )
  {
    case l1:
      x_->norm1( norm_view );
      break;
    case l2:
      x_->norm2( norm_view );
      break;
    default:
      x_->normInf( norm_view );
      break;
  }

  return norms[0];
}

//-----------------------------------------------------------------------------

auto Vector::min() const -> real
{
  dolfin_assert( not x_.is_null() );

  Teuchos::ArrayRCP< real const > arr = x_->getData( 0 );
  real min = *std::min_element( arr.get(), arr.get() + arr.size() );

  MPI::all_reduce_in_place< MPI::min >( min );

  return min;
}

//-----------------------------------------------------------------------------

auto Vector::max() const -> real
{
  dolfin_assert( not x_.is_null() );

  Teuchos::ArrayRCP< real const > arr = x_->getData( 0 );
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
  dolfin_assert( not x_.is_null() );

  x_->scale( a );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator/=( const real a ) -> Vector &
{
  dolfin_assert( not x_.is_null() );
  dolfin_assert( a != 0.0 );

  real const b = 1.0 / a;
  ( *this ) *= b;

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator*=( const GenericVector & y ) -> Vector &
{
  dolfin_assert( not x_.is_null() );

  trilinos::Vector const & Y_ = y.down_cast< trilinos::Vector >();
  dolfin_assert( not Y_.x_ghosted_.is_null() );

  x_->elementWiseMultiply( 1.0, *( x_->getVector( 0 ) ), *( Y_.x_ ), 0.0 );

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
    dolfin_assert( not v.x_.is_null() );
    dolfin_assert( not x_.is_null() );

    x_->assign( *v.x_ );
  }

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::operator=( real a ) -> Vector &
{
  dolfin_assert( not x_.is_null() );

  x_->putScalar( a );

  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::vec() const -> TPVectorPtr
{
  return x_;
}

//-----------------------------------------------------------------------------

auto Vector::factory() const -> LinearAlgebraFactory &
{
  return trilinos::Factory::instance();
}

//-----------------------------------------------------------------------------

auto Vector::clear() -> void
{
  if ( not x_.is_null() )
  {
    x_         = Teuchos::null;
    x_ghosted_ = Teuchos::null;
  }
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
