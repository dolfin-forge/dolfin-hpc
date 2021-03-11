// Copyright (C) 2004-2007 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScVector.h>

#include <dolfin/la/petsc/PETScFactory.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/PE.h>
#include <dolfin/math/basic.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
PETScVector::PETScVector()
  : Variable( "x", "a sparse vector" )

{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector( size_t N, bool distributed )
  : Variable( "x", "a sparse vector" )
  , x_( nullptr )
  , is_distributed_( false )
  , is_ghosted_( false )
{
  // Create PETSc vector
  init( N, distributed );
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector( Vec x )
  : Variable( "x", "a vector" )
  , x_( x )
  , is_distributed_( false )
  , is_ghosted_( false )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector( PETScVector const & v )
  : Variable( "x", "a vector" )
  , x_( nullptr )
  , is_distributed_( false )
  , is_ghosted_( false )
{
  *this = v;
}
//-----------------------------------------------------------------------------
PETScVector::~PETScVector()
{
  clear();
}
//-----------------------------------------------------------------------------
void PETScVector::init( size_t N, bool distributed )
{
  // Do not reallocate
  if ( x_ && this->local_size() == N )
  {
    return;
  }

  clear();

  PetscInt N_ = N;

  // Create vector
  if ( distributed && PE::size() > 1 )
  {
    is_distributed_ = true;
#ifdef DOLFIN_HAVE_MPI
    VecCreateMPI( MPI::DOLFIN_COMM, N_, PETSC_DETERMINE, &x_ );
#endif
  }
  else
  {
    VecCreate( PETSC_COMM_SELF, &x_ );
    VecSetSizes( x_, PETSC_DECIDE, N_ );
    VecSetFromOptions( x_ );
  }
  // Set all entries to zero
  PetscScalar a = 0.0;
  VecSet( x_, a );
}
//-----------------------------------------------------------------------------
void PETScVector::get( real * values ) const
{
#if PETSC_VERSION_MAJOR > 2
  dolfin_assert( x_ );
  real const * data = nullptr;
  VecGetArrayRead( x_, &data );
  dolfin_assert( data );
  PetscInt n;
  VecGetLocalSize( x_, &n );
  std::copy( data, data + n, values );
  VecRestoreArrayRead( x_, &data );
  dolfin_assert( x_ );
#else
  dolfin_assert( x_ );

  real * data = 0;
  VecGetArray( x_, &data );
  dolfin_assert( data );

  for ( size_t i = 0; i < local_size(); i++ )
    values[i] = data[i];
  VecRestoreArray( x_, &data );

  dolfin_assert( x_ );
#endif
}
//-----------------------------------------------------------------------------
void PETScVector::set( real * values )
{
  dolfin_assert( x_ );
  real * data = nullptr;
  VecGetArray( x_, &data );
  dolfin_assert( data );
  PetscInt n;
  VecGetLocalSize( x_, &n );
  std::copy( values, values + n, data );
  VecRestoreArray( x_, &data );
  dolfin_assert( x_ );
}
//-----------------------------------------------------------------------------
void PETScVector::add( real * values )
{
  dolfin_assert( x_ );
  PetscInt n;
  VecGetLocalSize( x_, &n );
  int * rows = new int[n];
  for ( int i = 0; i < n; i++ )
  {
    rows[i] = i;
  }
  VecSetValues( x_, n, rows, values, ADD_VALUES );
  delete[] rows;
}
//-----------------------------------------------------------------------------
void PETScVector::get( real * block, size_t m, const size_t * rows ) const
{
  dolfin_assert( x_ );

  PetscInt M_ = m;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );

  if ( is_ghosted_ )
  {
    PetscInt low, high;
    Vec xl;
    VecGetOwnershipRange( x_, &low, &high );
    VecGhostGetLocalForm( x_, &xl );

    std::vector< PetscInt > tmp( m );
    for ( size_t i = 0; i < m; i++ )
    {
      if ( rows_[i] < high && rows_[i] >= low )
      {
        tmp[i] = rows_[i] - low;
      }
      else
      {
        dolfin_assert( mapping_.size() > 0 );
        GhostMapping::const_iterator it = mapping_.find( rows[i] );
        dolfin_assert( mapping_.count( rows[i] ) > 0 );
        tmp[i] = it->second;
      }
    }
    VecGetValues( xl, M_, tmp.data(), block );
    VecGhostRestoreLocalForm( x_, &xl );
  }
  else
  {
    VecGetValues( x_, M_, rows_.data(), block );
  }
}
//-----------------------------------------------------------------------------
auto PETScVector::norm( VectorNormType type ) const -> real
{
  dolfin_assert( x_ );

  real value = 0.0;

  switch ( type )
  {
    case l1:
      VecNorm( x_, NORM_1, &value );
      break;
    case l2:
      VecNorm( x_, NORM_2, &value );
      break;
    default:
      VecNorm( x_, NORM_INFINITY, &value );
      break;
  }

  return value;
}
//-----------------------------------------------------------------------------
void PETScVector::pointwise( const GenericVector & x,
                             VectorPointwiseOp     op ) const
{

  const PETScVector & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  switch ( op )
  {
    case pw_min:
      VecPointwiseMin( x_, x_, v.x_ );
      break;
    case pw_max:
      VecPointwiseMax( x_, x_, v.x_ );
      break;
    case pw_mult:
      VecPointwiseMult( x_, x_, v.x_ );
      break;
    case pw_div:
      VecPointwiseDivide( x_, x_, v.x_ );
      break;
    case pw_maxabs:
      VecPointwiseMaxAbs( x_, x_, v.x_ );
      break;
    default:
      error( "Unknown operator" );
  }
}
//-----------------------------------------------------------------------------
void PETScVector::disp( size_t ) const
{
  section( "PETScVector" );
  if ( PE::size() > 1 && is_distributed_ )
  {
    VecView( x_, PETSC_VIEWER_STDOUT_WORLD );
  }
  else
  {
    VecView( x_, PETSC_VIEWER_STDOUT_SELF );
  }
  end();
}
//-----------------------------------------------------------------------------
void PETScVector::init_ghosted( size_t,
                                _ordered_set< size_t > &         indices,
                                _ordered_map< size_t, size_t > & map )
{
  if ( !is_distributed_ )
  {
    return;
  }

  if ( is_ghosted_ )
  {
    apply();
  }

  PetscInt local_size, size, low, high;
  VecGetSize( x_, &size );
  VecGetLocalSize( x_, &local_size );
  VecGetOwnershipRange( x_, &low, &high );

  mapping_.clear();

  std::vector < PetscInt > rows( local_size );
  std::vector <real >      values( local_size );
  for ( int i = 0; i < local_size; i++ )
  {
    rows[i]           = low + i;
    mapping_[low + i] = i;
  }

  VecGetValues( x_, local_size, rows.data(), values.data() );

  if ( is_ghosted_ && map.size() > 0 )
  {
    // dolfin_assert(map.size() > 0);
    for ( int i = 0; i < local_size; i++ )
    {
      // dolfin_assert(map.count(low + i) > 0);
      rows[i] = map[low + i];
    }
  }

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  VecDestroy( &x_ );
#else
  VecDestroy( x_ );
#endif

  std::vector< PetscInt > ghost_indices;
  PetscInt                num_ghost = local_size;

  for ( size_t const & index : indices )
  {
    if ( index < static_cast< size_t >( low )
         or index >= static_cast< size_t >( high ) )
    {
      ghost_indices.push_back( index );
      mapping_[index] = num_ghost++;
    }
  }

#ifdef DOLFIN_HAVE_MPI
  num_ghost = ghost_indices.size();
  VecCreateGhost( MPI::DOLFIN_COMM, local_size, size,
                  num_ghost, ghost_indices.data(), &x_ );
#endif
  VecSetValues( x_, local_size, rows.data(), values.data(), INSERT_VALUES );

  is_ghosted_ = true;
  apply();
}

//-----------------------------------------------------------------------------
auto PETScVector::factory() const -> LinearAlgebraFactory &
{
  return PETScFactory::instance();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */
