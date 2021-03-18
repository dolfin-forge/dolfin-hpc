// Copyright (C) 2020 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosVector.h>

#include <Teuchos_OrdinalTraits.hpp>

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

Vector::Vector( TPVector x )
  : Variable( "x", "a vector" )
{
  // Do nothing
}

//-----------------------------------------------------------------------------

Vector::Vector( Vector const & v )
  : Variable( "x", "a vector" )
  // , x_( Teuchos::rcp( new TPVector::element_type( *v.x_ ) ) )
  // , contigMap( Teuchos::rcp( new TPMap::element_type( *v.contigMap ) ) )
{
  *this = v;
}

//-----------------------------------------------------------------------------

Vector::~Vector()
{
  clear();
}

//-----------------------------------------------------------------------------

Vector * Vector::copy() const
{
  return new Vector( *this );
}

//-----------------------------------------------------------------------------

void Vector::zero()
{
  // delegate to assignment operator with 0.0
  *this = 0.0;
}

//-----------------------------------------------------------------------------
void Vector::apply( FinalizeType finaltype )
{
  // FIXME
}

//-----------------------------------------------------------------------------
void Vector::disp( size_t precision ) const
{
  // FIXME
}

//-----------------------------------------------------------------------------

size_t Vector::size() const
{
  // if ( not contigMap.is_null() )
  // {
  //   return static_cast< size_t >( contigMap->getGlobalNumElements() );
  // }

  return 0;
}

//-----------------------------------------------------------------------------

size_t Vector::local_size() const
{
  // if ( not contigMap.is_null() )
  // {
  //   return static_cast< size_t >( contigMap->getNodeNumElements() );
  // }

  return 0;
}

//-----------------------------------------------------------------------------

size_t Vector::offset() const
{
  // dolfin_assert( not x_.is_null() );
  // // this is wrong!!!
  // int low = contigMap->getMinLocalIndex();
  // return static_cast< size_t >( low );
}

//-----------------------------------------------------------------------------

void Vector::init( size_t N )
{
  init( N, true );
}

//-----------------------------------------------------------------------------

void Vector::init( size_t N, bool distributed )
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

  size_t const rank = MPI::rank();
  size_t const size = MPI::size();

  // Compute number of items per process and remainder
  size_t const n = N / size;
  size_t const r = N % size;

  // Mapping across processes
  size_t Nlocal = ( rank < r )
                  ? ( rank * ( n + 1 ) + n + 1 ) - ( rank * ( n + 1 ) )
                  : ( rank * n + r + n ) - ( rank * n + r );


#if DEBUG
  // check if we get N elements across all processes
  {
    size_t Ntest = 0;
    MPI::all_reduce< MPI::sum >( Nlocal, Ntest );
    dolfin_assert( N == Ntest );
  }
#endif

  Teuchos::RCP< const Teuchos::MpiComm< int > > _comm(
    new Teuchos::MpiComm< int >( Teuchos::MpiComm< int >( MPI::DOLFIN_COMM ) ) );

  // Teuchos::RCP< TPMap > _map( new TPMap( N, Nlocal, 0, _comm ) );
  // Teuchos::RCP< TPMap > _ghost_map;

  // std::vector< int > local_to_global_map;

  // // Save a map for the ghosting of values on other processes
  // if ( local_to_global_map.size() != 0 )
  // {
  //   const Teuchos::ArrayView< const int > local_indices( local_to_global_map );
  //   _ghost_map = Teuchos::rcp( new TPMap( N, local_indices, 0, _comm ) );
  // }
  // else
  // {
  //   _ghost_map = _map;
  // }

  // // Vector - create with overlap
  // x_ghosted_ = Teuchos::rcp( new TPVector( _ghost_map, 1 ) );

  // // make sure we actually got a non-empty vector
  // dolfin_assert( not x_ghosted_.is_null() );

  // // Get a modifiable view into the ghosted vector
  // x_ = x_ghosted_->offsetViewNonConst( _map, 0 );
}

//-----------------------------------------------------------------------------

void Vector::init_ghosted( size_t                           n,
                           _ordered_set< size_t > &         indices,
                           _ordered_map< size_t, size_t > & map )
{
  // FIXME
}

//-----------------------------------------------------------------------------

void Vector::get( real * values ) const
{
  // #if PETSC_VERSION_MAJOR > 2
  //   dolfin_assert(x_);
  //   real const* data = nullptr;
  //   VecGetArrayRead(x_, &data);
  //   dolfin_assert(data);
  //   PetscInt n;
  //   VecGetLocalSize(x_, &n);
  //   std::copy(data, data + n, values);
  //   VecRestoreArrayRead(x_, &data);
  //   dolfin_assert(x_);
  // #else
  //   dolfin_assert(x_);

  //   real* data = 0;
  //   VecGetArray(x_, &data);
  //   dolfin_assert(data);

  //   for (size_t i = 0; i < local_size(); i++)
  //     values[i] = data[i];
  //   VecRestoreArray(x_, &data);

  //   dolfin_assert(x_);
  // #endif
}

//-----------------------------------------------------------------------------

void Vector::set( real * values )
{
  // dolfin_assert(x_);
  // real* data = nullptr;
  // VecGetArray(x_, &data);
  // dolfin_assert(data);
  // PetscInt n;
  // VecGetLocalSize(x_, &n);
  // std::copy(values, values + n, data);
  // VecRestoreArray(x_, &data);
  // dolfin_assert(x_);
}

//-----------------------------------------------------------------------------

void Vector::add( real * values )
{
  // dolfin_assert(x_);
  // PetscInt n;
  // VecGetLocalSize(x_, &n);
  // int * rows = new int[n];
  // for (int i = 0; i < n; i++) { rows[i] = i; }
  // VecSetValues(x_, n, rows, values, ADD_VALUES);
  // delete[] rows;
}

//-----------------------------------------------------------------------------

void Vector::get( real * block, size_t m, const size_t * rows ) const
{
  // dolfin_assert(x_);

  // if (is_ghosted_)
  // {
  //   int low, high;
  //   Vec xl;
  //   VecGetOwnershipRange(x_, &low, &high);
  //   VecGhostGetLocalForm(x_, &xl);

  //   int *tmp = new int[m];
  //   for (size_t i = 0; i < m; i++)
  //   {
  //     if ((int) rows[i] < high && (int) rows[i] >= low)
  //     {
  //       tmp[i] = rows[i] - low;
  //     }
  //     else
  //     {
  //       dolfin_assert(mapping_.size() > 0);
  //       GhostMapping::const_iterator it = mapping_.find(rows[i]);
  //       dolfin_assert(mapping_.count(rows[i]) > 0);
  //       tmp[i] = it->second;
  //     }
  //   }
  //   VecGetValues(xl, static_cast<int>(m), tmp, block);
  //   VecGhostRestoreLocalForm(x_, &xl);

  //   delete[] tmp;
  // }
  // else
  // {
  //   VecGetValues(x_, static_cast<int>(m),
  //                reinterpret_cast<int*>(const_cast<size_t*>(rows)), block);
  // }
}

//-----------------------------------------------------------------------------

void Vector::set( const real * block, size_t m, const size_t * rows )
{
  // dolfin_assert( x_ );

  // const TPVector::element_type::global_ordinal_type gblRow = 2;
  // A(0, 0:1) = [2, -1]
  // if ( gblRow == 0 )
  // {
  //   x_->insertGlobalValues( gblRow, gblRow, 1.0 );
  // }

  // VecSetValues( x_,
  //               static_cast< int >( m ),
  //               reinterpret_cast< int * >( const_cast< size_t * >( rows ) ),
  //               block,
  //               INSERT_VALUES );
}

//-----------------------------------------------------------------------------

void Vector::add( const real * block, size_t m, const size_t * rows )
{
  // dolfin_assert( x_ );
  // VecSetValues(x_, static_cast<int>(m),
  //              reinterpret_cast<int*>(const_cast<size_t*>(rows)), block,
  //              ADD_VALUES);
}

//-----------------------------------------------------------------------------

void Vector::axpy( real a, const GenericVector & y )
{
  // dolfin_assert( x_.is_null() );

  // Vector const & v = y.down_cast< Vector >();
  // dolfin_assert( v.x_.is_null() );

  // if ( size() != v.size() )
  // {
  //   error( "The vectors must be of the same size to apply AXPY." );
  // }

  // // x_ = a * x_ + v.x_ * 1.0
  // x_->update( a, *v.vec(), 1.0 );
}

//-----------------------------------------------------------------------------

void Vector::axpby( real a, const GenericVector & x, real b )
{
  // FIXME
}

//-----------------------------------------------------------------------------

void Vector::waxpy( real a, const GenericVector & x, const GenericVector & y )
{
  // FIXME
}

//-----------------------------------------------------------------------------

void Vector::axpbypcz( real                  a,
                       const GenericVector & x,
                       real                  b,
                       const GenericVector & y,
                       real                  c )
{
  // FIXME
}

//-----------------------------------------------------------------------------

real Vector::inner( const GenericVector & y ) const
{
  // dolfin_assert( x_.is_null() );

  // Vector const & v = y.down_cast< Vector >();
  // dolfin_assert( v.x_.is_null() );

  // return x_->dot( *v.vec() );
}

//-----------------------------------------------------------------------------

real Vector::norm( VectorNormType type ) const
{
  // dolfin_assert( not x_.is_null() );

  real value = 0.0;

  // switch ( type )
  // {
  //   case l1:
  //     value = x_->norm1();
  //     break;
  //   case l2:
  //     value = x_->norm2();
  //     break;
  //   default:
  //     value = x_->normInf();
  //     break;
  // }

  return value;
}

//-----------------------------------------------------------------------------

real Vector::min() const
{
  // dolfin_assert( not x_.is_null() );

  real value = 0.0;

  // TODO

  return value;
}

//-----------------------------------------------------------------------------

real Vector::max() const
{
  // dolfin_assert( not x_.is_null() );

  real value = 0.0;

  // TODO

  return value;
}

//-----------------------------------------------------------------------------

void Vector::pointwise( const GenericVector & x,
                        VectorPointwiseOp     op ) const
{
}

//-----------------------------------------------------------------------------

Vector & Vector::operator*=( const real a )
{
  // dolfin_assert( x_.is_null() );
  // x_->scale( a );
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator/=( const real a )
{
  // dolfin_assert( x_.is_null() );
  dolfin_assert( a != 0.0 );
  // x_->scale( 1.0 / a );
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator*=( const GenericVector & y )
{
  // dolfin_assert( x_.is_null() );
  // Vector const & v = y.down_cast< Vector >();
  // dolfin_assert( v.x_.is_null() );

  // if ( size() != v.size() )
  // {
  //   error(
  //     "Vectors must have the same size for componentwise multiplication." );
  // }

  // TODO

  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator+=( const GenericVector & x )
{
  // this->axpy( 1.0, x );
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator-=( const GenericVector & x )
{
  // this->axpy( -1.0, x );
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator=( const GenericVector & v )
{
  // *this = v.down_cast< Vector >();
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator=( Vector const & v )
{
  // if ( &v != this )
  // {
  //   dolfin_assert( not v.x_.is_null() );
  //   init( v.local_size(), v.x_->isDistributed() );
  //   x_->assign( *v.x_ );
  // }
  return *this;
}

//-----------------------------------------------------------------------------

Vector & Vector::operator=( real a )
{
  // dolfin_assert( not x_.is_null() );
  // x_->putScalar( a );
  return *this;
}

//-----------------------------------------------------------------------------

auto Vector::vec() const -> TPVectorPtr
{
  return x_;
}

//-----------------------------------------------------------------------------

LinearAlgebraFactory & Vector::factory() const
{
  // FIXME
}

//-----------------------------------------------------------------------------

void Vector::clear()
{
  // if ( not x_.is_null() )
  // {
  //   x_        = Teuchos::null;
  //   contigMap = Teuchos::null;
  // }
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
