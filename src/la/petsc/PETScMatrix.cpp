// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScMatrix.h>

#include <dolfin/common/system.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/petsc/PETScFactory.h>
#include <dolfin/la/petsc/PETScVector.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <iomanip>
#include <iostream>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix()
  : Variable( "A", "a sparse matrix" )
{
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix( Mat A )
  : Variable( "A", "a sparse matrix" )
  , A( A )
  , AA_sub( nullptr )
  , is_distributed_( false )
  , rstart_( 0 )
  , rend_( 0 )
{
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix( size_t M, size_t N, bool distributed )
  : Variable( "A", "a sparse matrix" )
  , A( nullptr )
  , AA_sub( nullptr )
  , is_distributed_( false )
  , rstart_( 0 )
  , rend_( 0 )
{
  init( M, N, distributed );
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix( const PETScMatrix & A )
  : Variable( "A", "PETSc matrix" )
  , A( nullptr )
  , AA_sub( nullptr )
  , is_distributed_( false )
  , rstart_( 0 )
  , rend_( 0 )
{
  *this = A;
}
//-----------------------------------------------------------------------------
PETScMatrix::~PETScMatrix()
{
  clear();
}
//-----------------------------------------------------------------------------
void PETScMatrix::clear()
{
  if ( A )
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    MatDestroy( &A );
#else
    MatDestroy( A );
#endif
  }

  if ( AA_sub )
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    MatDestroyMatrices( 1, &AA_sub );
#else
    MatDestroyMatrices( 1, &AA_sub );
#endif
    AA_sub = nullptr;
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::init( size_t M, size_t N, bool distributed )
{
  // Free previously allocated memory if necessary
  clear();

  PetscInt M_ = M;
  PetscInt N_ = N;

  // Create a sparse matrix in compressed row format
  if ( dolfin::MPI::size() > 1 && distributed )
  {
    is_distributed_ = true;
    // Create PETSc parallel matrix with a guess for number of diagonal and
    /// number of off-diagonal non-zeroes.
    // Note that guessing too high leads to excessive memory usage.
#ifdef DOLFIN_HAVE_MPI
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3
    MatCreateAIJ( dolfin::MPI::DOLFIN_COMM, M_, N_,
                  PETSC_DETERMINE, PETSC_DETERMINE,
                  120, PETSC_NULL,
                  120, PETSC_NULL,
                  &A );
#else
    MatCreateMPIAIJ( dolfin::MPI::DOLFIN_COMM, M_, N_,
                     PETSC_DETERMINE, PETSC_DETERMINE,
                     120, PETSC_NULL,
                     120, PETSC_NULL,
                     &A );
#endif
#endif
  }
  else
  {
    // Create PETSc sequential matrix with a guess for number of non-zeroes
    MatCreateSeqAIJ( PETSC_COMM_SELF, M_, N_, 50, PETSC_NULL, &A );

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
    MatSetOption( A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE );
#else
    MatSetOption( A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE );
#endif
#else
    MatSetOption( A, MAT_KEEP_ZEROED_ROWS );
#endif
    MatSetFromOptions( A );
  }

  MatGetOwnershipRange( A, &rstart_, &rend_ );
}
//-----------------------------------------------------------------------------
void PETScMatrix::init( size_t M, size_t N, std::vector< size_t > const & nz )
{
  dolfin_assert( is_distributed_ == false );

  // Free previously allocated memory if necessary
  clear();

  PetscInt M_ = M;
  PetscInt N_ = N;

  std::vector< PetscInt > nz_( nz.begin(), nz.end() );

  // Create a sparse matrix in compressed row format
  MatCreate( PETSC_COMM_SELF, &A );
  MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, M_, N_ );
  MatSetType( A, MATSEQAIJ );

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
  MatSetOption( A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE );
#else
  MatSetOption( A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE );
#endif
#else
  MatSetOption( A, MAT_KEEP_ZEROED_ROWS );
#endif
  MatSetFromOptions( A );
  // MatSeqAIJSetPreallocation( A, PETSC_DEFAULT, ( int * ) nz );
  MatSeqAIJSetPreallocation( A, PETSC_DEFAULT, nz_.data() );
  MatZeroEntries( A );

  MatGetOwnershipRange( A, &rstart_, &rend_ );
}
//-----------------------------------------------------------------------------
void PETScMatrix::init( size_t M, size_t N,
                        std::vector< size_t > const & d_nzrow,
                        std::vector< size_t > const & o_nzrow )
{
  dolfin_assert( is_distributed_ == true );

  // Free previously allocated memory if necessary
  clear();

  // Create PETSc parallel matrix with a guess for number of diagonal (50 in
  // this case) and number of off-diagonal non-zeroes (50 in this case). Note
  // that guessing too high leads to excessive memory usage. In order to not
  // waste any memory one would need to specify d_nnz and o_nnz.
#ifdef DOLFIN_HAVE_MPI

  PetscInt M_ = M;
  PetscInt N_ = N;

  std::vector< PetscInt > d_nzrow_( d_nzrow.begin(), d_nzrow.end() );
  std::vector< PetscInt > o_nzrow_( o_nzrow.begin(), o_nzrow.end() );

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3
  MatCreateAIJ( MPI::DOLFIN_COMM, M_, N_,
                PETSC_DETERMINE, PETSC_DETERMINE,
                PETSC_DETERMINE, d_nzrow_.data(),
                PETSC_DETERMINE, o_nzrow_.data(),
                &A );
#else
  MatCreateMPIAIJ( MPI::DOLFIN_COMM, M_, N_,
                   PETSC_DETERMINE, PETSC_DETERMINE,
                   PETSC_DETERMINE, d_nzrow_.data(),
                   PETSC_DETERMINE, o_nzrow_.data(),
                   &A );
#endif
#else
  MAYBE_UNUSED( M );
  MAYBE_UNUSED( N );
  MAYBE_UNUSED( d_nzrow );
  MAYBE_UNUSED( o_nzrow );
#endif

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
  MatSetOption( A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE );
#else
  MatSetOption( A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE );
#endif
#else
  MatSetOption( A, MAT_KEEP_ZEROED_ROWS );
#endif
  MatSetFromOptions( A );
  MatZeroEntries( A );

  MatGetOwnershipRange( A, &rstart_, &rend_ );
}
//-----------------------------------------------------------------------------
void PETScMatrix::init( const GenericSparsityPattern & sparsity_pattern )
{
  if ( sparsity_pattern.is_distributed() )
  {
    SparsityPattern const & spattern =
      reinterpret_cast< SparsityPattern const & >( sparsity_pattern );

    is_distributed_        = true;
    size_t p               = dolfin::MPI::rank();
    size_t local_size      = spattern.range_size( p );
    size_t spattern_size_0 = spattern.size( 0 );
    size_t spattern_size_1 = spattern.size( 1 );

    std::vector< size_t > d_nzrow( local_size );
    std::vector< size_t > o_nzrow( local_size );
    spattern.numNonZeroPerRow( p, d_nzrow.data(), o_nzrow.data() );

    const_cast< SparsityPattern & >( spattern ).clear();
    init( spattern_size_0, spattern_size_1, d_nzrow, o_nzrow );
  }
  else
  {
    SparsityPattern const & spattern =
      reinterpret_cast< SparsityPattern const & >( sparsity_pattern );

    std::vector< size_t > nzrow( spattern.size( 0 ) );
    spattern.numNonZeroPerRow( nzrow.data() );

    init( spattern.size( 0 ), spattern.size( 1 ), nzrow );
  }
}
//-----------------------------------------------------------------------------
auto PETScMatrix::copy() const -> PETScMatrix *
{
  PETScMatrix * mcopy = new PETScMatrix();
  MatDuplicate( A, MAT_COPY_VALUES, &( mcopy->A ) );

  return mcopy;
}
//-----------------------------------------------------------------------------
auto PETScMatrix::norm( std::string norm_type ) const -> real
{
  real norm;
  if ( norm_type == "l1" )
  {
    MatNorm( A, NORM_1, &norm );
  }
  else if ( norm_type == "linf" )
  {
    MatNorm( A, NORM_INFINITY, &norm );
  }
  else if ( norm_type == "frobenius" )
  {
    MatNorm( A, NORM_FROBENIUS, &norm );
  }
  else
  {
    error( "Unknown norm type in uPETScMatrix." );
    return 0.0;
  }
  return norm;
}
//-----------------------------------------------------------------------------
void PETScMatrix::getrow( size_t                  row,
                          std::vector< size_t > & columns,
                          std::vector< real > &   values ) const
{
  const PetscInt *    cols  = nullptr;
  const PetscScalar * vals  = nullptr;
  PetscInt            ncols = 0;
  if ( row >= static_cast< size_t >( rstart_ )
       && row < static_cast< size_t >( rend_ ) )
  {
    // Local row
    MatGetRow( A, row, &ncols, &cols, &vals );
    columns.assign( cols, cols + ncols );
    values.assign( vals, vals + ncols );
    MatRestoreRow( A, row, &ncols, &cols, &vals );
  }
  else
  {
    // Non-local row
    if ( AA_sub )
    {
      error( "Trying to get ghosted row %d but no SubMatrices defined.\n"
             "Ownership range %d, %d", row, rstart_, rend_ );
    }

    _map< int, int >::const_iterator it = mapping_.find( row );
    MatGetRow( AA_sub[0], it->second, &ncols, &cols, &vals );
    columns.assign( cols, cols + ncols );
    values.assign( vals, vals + ncols );
    MatRestoreRow( AA_sub[0], it->second, &ncols, &cols, &vals );
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::getrows_offproc( _ordered_set< size_t > const & rows )
{
  if ( !is_distributed_ )
  {
    return;
  }

  int * _cols = new int[size( 0 )];
  int * _rows = new int[rows.size()];

  mapping_.clear();

  PetscInt i = 0;
  for ( size_t const & row : rows )
  {
    if ( row < static_cast< size_t >( rstart_ )
         && row >= static_cast< size_t >( rend_ ) )
    {
      _rows[i]      = row;
      mapping_[row] = i++;
    }
  }

  for ( size_t j = 0; j < size( 0 ); ++j )
  {
    _cols[j] = j;
  }

  IS irow, icol;
  PetscInt j = size( 0 );

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral( PETSC_COMM_SELF, i, _rows, PETSC_COPY_VALUES, &irow );
  ISCreateGeneral( PETSC_COMM_SELF, j, _cols, PETSC_COPY_VALUES, &icol );

  if ( !AA_sub )
  {
    MatCreateSubMatrices( A, 1, &irow, &icol, MAT_INITIAL_MATRIX, &AA_sub );
  }
  else
  {
    MatCreateSubMatrices( A, 1, &irow, &icol, MAT_REUSE_MATRIX, &AA_sub );
  }

  ISDestroy( &irow );
  ISDestroy( &icol );
#else
  ISCreateGeneral( PETSC_COMM_SELF, i, _rows, &irow );
  ISCreateGeneral( PETSC_COMM_SELF, j, _cols, &icol );

  if ( !AA_sub )
  {
    MatGetSubMatrices( A, 1, &irow, &icol, MAT_INITIAL_MATRIX, &AA_sub );
  }
  else
  {
    MatGetSubMatrices( A, 1, &irow, &icol, MAT_REUSE_MATRIX, &AA_sub );
  }

  ISDestroy( irow );
  ISDestroy( icol );
#endif

  delete[] _cols;
  delete[] _rows;
}
//-----------------------------------------------------------------------------
void PETScMatrix::zero( size_t m, size_t const * rows )
{
  IS is            = nullptr;
  PetscScalar null = 0.0;
  PetscInt m_      = m;
  std::vector< PetscInt > rows_( rows, rows + m );

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral( PETSC_COMM_SELF, m_, rows_.data(), PETSC_COPY_VALUES, &is );
  MatZeroRowsIS( A, is, null, PETSC_NULL, PETSC_NULL );
  ISDestroy( &is );
#else
  ISCreateGeneral( PETSC_COMM_SELF, m_, rows_.data(), &is );
  MatZeroRowsIS( A, is, null );
  ISDestroy( is );
#endif
}
//-----------------------------------------------------------------------------
void PETScMatrix::ident( size_t m, size_t const * rows )
{
  IS is           = nullptr;
  PetscScalar one = 1.0;
  PetscInt m_     = m;
  std::vector< PetscInt > rows_( rows, rows + m );

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral( PETSC_COMM_SELF, m_, rows_.data(), PETSC_COPY_VALUES, &is );
  MatZeroRowsIS( A, is, one, PETSC_NULL, PETSC_NULL );
  ISDestroy( &is );
#else
  ISCreateGeneral( PETSC_COMM_SELF, m_, rows_.data(), &is );
  MatZeroRowsIS( A, is, one );
  ISDestroy( is );
#endif
}
//-----------------------------------------------------------------------------
void PETScMatrix::mult( const GenericVector & x,
                        GenericVector &       y,
                        bool                  transposed ) const
{
  const PETScVector & xx = x.down_cast< PETScVector >();
  PETScVector &       yy = y.down_cast< PETScVector >();

  if ( transposed )
  {
    if ( size( 0 ) != xx.size() )
    {
      error(
        "Matrix and vector dimensions mismatch for matrix-vector product." );
    }
    yy.init( xx.local_size() );
    MatMultTranspose( A, xx.vec(), yy.vec() );
  }
  else
  {
    if ( size( 1 ) != xx.size() )
    {
      error(
        "Matrix and vector dimensions mismatch for matrix-vector product." );
    }
    yy.init( xx.local_size() );
    MatMult( A, xx.vec(), yy.vec() );
  }
}
//-----------------------------------------------------------------------------
auto PETScMatrix::norm( const Norm type ) const -> real
{
  real value = 0.0;
  switch ( type )
  {
    case l1:
      MatNorm( A, NORM_1, &value );
      break;
    case linf:
      MatNorm( A, NORM_INFINITY, &value );
      break;
    case frobenius:
      MatNorm( A, NORM_FROBENIUS, &value );
      break;
    default:
      error( "Unknown norm type." );
      break;
  }
  return value;
}
//-----------------------------------------------------------------------------
void PETScMatrix::disp( size_t ) const
{
  section( "PETScMatrix" );
  MatInfo info;

  section( "Local" );
  MatGetInfo( A, MAT_LOCAL, &info );
  print( info );
  end();
  skip();

  if ( is_distributed_ )
  {
    section( "Global" );
    MatGetInfo( A, MAT_GLOBAL_SUM, &info );
    print( info );
    end();
    skip();
  }

  if ( verbose() == 0 )
  {
    return;
  }

  if ( is_distributed_ )
  {

    MatView( A, PETSC_VIEWER_STDOUT_WORLD );
  }
  else
  {
    PetscViewerPushFormat( PETSC_VIEWER_STDOUT_SELF,
                           PETSC_VIEWER_ASCII_MATLAB );
#ifdef BLOCKED
    PetscObjectSetName( ( PetscObject ) A, " Ab" );
#else
    PetscObjectSetName( ( PetscObject ) A, " Anb" );
#endif
    MatView( A, PETSC_VIEWER_STDOUT_SELF );
  }
  end();
}
//-----------------------------------------------------------------------------
auto PETScMatrix::factory() const -> LinearAlgebraFactory &
{
  return PETScFactory::instance();
}
//-----------------------------------------------------------------------------
void PETScMatrix::print( MatInfo const & info ) const
{
  message( "%24s : %lu", "block_size", uidx( info.block_size ) );
  message( "%24s : %lu", "nz_allocated", uidx( info.nz_allocated ) );
  message( "%24s : %lu", "nz_used", uidx( info.nz_used ) );
  message( "%24s : %lu", "nz_unneeded", uidx( info.nz_unneeded ) );
  message( "%24s : %s", "memory", human_readable( info.memory ).c_str() );
  message( "%24s : %lu", "assemblies", uidx( info.assemblies ) );
  message( "%24s : %lu", "mallocs", uidx( info.mallocs ) );
  message( "%24s : %lu", "fill_ratio_given", uidx( info.fill_ratio_given ) );
  message( "%24s : %lu", "fill_ratio_needed", uidx( info.fill_ratio_needed ) );
  message( "%24s : %lu", "factor_mallocs", uidx( info.factor_mallocs ) );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */
