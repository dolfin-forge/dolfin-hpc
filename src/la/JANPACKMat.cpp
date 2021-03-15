// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#include <string>

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/la/JANPACKFactory.h>
#include <dolfin/la/JANPACKMat.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/log/dolfin_log.h>

#ifdef HAVE_JANPACK

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#include <dolfin/main/MPI.h>
#endif

#include <janpack/spmv.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

JANPACKMat::JANPACKMat()
  : Variable( "A", "JANPACK matrix" )
#ifdef HAVE_JANPACK_MPI
  , A( &AA )
#endif
{
  // TODO: call JANPACK_Init or something?
}

//-----------------------------------------------------------------------------

JANPACKMat::JANPACKMat( size_t M, size_t N )
  : Variable( "A", "JANPACK matrix" )
{
  // TODO: call JANPACK_Init or something?
  // Create JANPACK matrix
  init( M, N );
}

//-----------------------------------------------------------------------------

JANPACKMat::JANPACKMat( const JANPACKMat & A )
  : Variable( "A", "JANPACK matrix" )
{
  error( "Not implemented" );
}

//-----------------------------------------------------------------------------

JANPACKMat::~JANPACKMat()
{
  jp_mat_free( A );
}

//-----------------------------------------------------------------------------

void JANPACKMat::init( size_t M, size_t N )
{
  jp_mat_init( A, M, N );

  // Allow for zeros to be inserted into the matrix
  jp_mat_setopt( A, JP_MAT_ZEROS );
}

//-----------------------------------------------------------------------------

void JANPACKMat::init( size_t M, size_t N, bool distributed )
{
  init( M, N );
}

//-----------------------------------------------------------------------------

void JANPACKMat::init( const GenericSparsityPattern & sparsity_pattern )
{

  const SparsityPattern & spattern =
    reinterpret_cast< const SparsityPattern & >( sparsity_pattern );

  init( spattern.size( 0 ), spattern.size( 1 ) );
}

//-----------------------------------------------------------------------------

JANPACKMat * JANPACKMat::copy() const
{

  error( "JANPACKMat::copy not yet implemented." );

  JANPACKMat * mcopy = new JANPACKMat();

  return mcopy;
}

//-----------------------------------------------------------------------------

dolfin::size_t JANPACKMat::size( size_t dim ) const
{
  dolfin_assert( A );
  size_t32_t M = 0;
  size_t32_t N = 0;
  jp_mat_size( const_cast< jp_mat_type * >( A ), &M, &N );
  return ( dim == 0 ? M : N );
}

//-----------------------------------------------------------------------------

dolfin::size_t JANPACKMat::nz() const
{
  return jp_mat_nz( const_cast< jp_mat_type * >( A ), 1 );
}

//-----------------------------------------------------------------------------

void JANPACKMat::get( real *         block,
                      size_t         m,
                      const size_t * rows,
                      size_t         n,
                      const size_t * cols ) const
{
  dolfin_assert( A );

  // Not yet implemented
  error( "JANPACKMat::get not yet implemented." );
}

//-----------------------------------------------------------------------------

void JANPACKMat::set( const real *   block,
                      size_t         m,
                      const size_t * rows,
                      size_t         n,
                      const size_t * cols )
{
  dolfin_assert( A );

  jp_mat_set_block( A,
                    m,
                    const_cast< size_t * >( rows ),
                    n,
                    const_cast< size_t * >( cols ),
                    const_cast< real * >( block ) );
}

//-----------------------------------------------------------------------------

void JANPACKMat::add( const real *   block,
                      size_t         m,
                      const size_t * rows,
                      size_t         n,
                      const size_t * cols )
{
  dolfin_assert( A );

  jp_mat_add_block( A,
                    m,
                    const_cast< size_t * >( rows ),
                    n,
                    const_cast< size_t * >( cols ),
                    const_cast< real * >( block ) );
}

//-----------------------------------------------------------------------------

real JANPACKMat::norm( std::string norm_type ) const
{
  error( "Not implemented." );
  return 0.0;
}

//-----------------------------------------------------------------------------

void JANPACKMat::zero()
{
  dolfin_assert( A );
  jp_mat_zero( A );
}

//-----------------------------------------------------------------------------

void JANPACKMat::apply( FinalizeType finaltype )
{

  jp_mat_finalize( A );
}

//-----------------------------------------------------------------------------

void JANPACKMat::disp( size_t precision ) const
{
#ifndef HAVE_JANPACK_MPI
  jp_mat_print( const_cast< jp_mat_type * >( A ) );
#endif
}

//-----------------------------------------------------------------------------

void JANPACKMat::ident( size_t m, const size_t * rows )
{
  dolfin_assert( A );

  for ( size_t i = 0; i < m; i++ )
  {
    // jp_mat_ident(A, rows[i]);
    jp_mat_zero_row( A, rows[i] );
    jp_mat_insert( A, rows[i], rows[i], 1.0 );
  }
}

//-----------------------------------------------------------------------------

void JANPACKMat::zero( size_t m, const size_t * rows )
{
  dolfin_assert( A );

  for ( size_t i = 0; i < m; i++ )
    jp_mat_zero_row( A, rows[i] );
}

//-----------------------------------------------------------------------------

void JANPACKMat::mult( const GenericVector & x,
                       GenericVector &       y,
                       bool                  transposed ) const
{
  dolfin_assert( A );
  const JANPACKVec & xx = x.down_cast< JANPACKVec >();
  JANPACKVec &       yy = y.down_cast< JANPACKVec >();
  if ( transposed )
    yy.init( xx.local_size() );
  else
    yy.init( xx.local_size() );

  jp_spmv( const_cast< jp_mat_type * >( A ), xx.vec(), yy.vec() );
}

//-----------------------------------------------------------------------------

void JANPACKMat::getrow( size_t                  row,
                         std::vector< size_t > & columns,
                         std::vector< real > &   values ) const
{

  size_t n, *c = 0;
  real * v = 0;

  c = new size_t[size( 0 )];
  v = new real[size( 0 )];

  jp_mat_getrow( const_cast< jp_mat_type * >( A ), row, c, v, &n );

  // Assign values to std::vectors
  columns.assign( reinterpret_cast< size_t * >( c ),
                  reinterpret_cast< size_t * >( c + n ) );
  values.assign( v, v + n );

  delete[] c;
  delete[] v;
}

//-----------------------------------------------------------------------------

void JANPACKMat::setrow( size_t                        row,
                         const std::vector< size_t > & columns,
                         const std::vector< real > &   values )
{
  set( &values[0], 1, &row, columns.size(), &columns[0] );
}

//-----------------------------------------------------------------------------

LinearAlgebraFactory & JANPACKMat::factory() const
{
  return JANPACKFactory::instance();
}

//-----------------------------------------------------------------------------

jp_mat_type * JANPACKMat::mat() const
{
  return const_cast< jp_mat_type * >( A );
}

//-----------------------------------------------------------------------------

const JANPACKMat & JANPACKMat::operator*=( real a )
{
  dolfin_assert( A );
  error( "Not implemented." );
  return *this;
}

//-----------------------------------------------------------------------------

const JANPACKMat & JANPACKMat::operator/=( real a )
{
  dolfin_assert( A );
  error( "Not implemented." );
  return *this;
}

//-----------------------------------------------------------------------------

const GenericMatrix & JANPACKMat::operator=( const GenericMatrix & A )
{
  //  jp_mat_copy(A.down_cast<JANPACKMat>().A, this->A);
  error( "Not implemented." );
  return *this;
}

//-----------------------------------------------------------------------------

void JANPACKMat::dup( const JANPACKMat & A )
{
  size_t range[2];
  size_t m;

  jp_mat_range( const_cast< jp_mat_type * >( A.mat() ), &range[0], &range[1] );
  m = range[1] - range[0];
  init( m, m );
}

//-----------------------------------------------------------------------------

#endif
