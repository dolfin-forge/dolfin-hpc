// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosMatrix.h>

#include <dolfin/la/trilinos/TrilinosFactory.h>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

Matrix::Matrix()
{
  // FIXME
}

//-----------------------------------------------------------------------------

Matrix::Matrix( size_t M, size_t N, bool distributed )
{
  // FIXME
}

//-----------------------------------------------------------------------------

Matrix::Matrix( const Matrix & A )
{
  // FIXME
}

//-----------------------------------------------------------------------------

Matrix::~Matrix()
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( const GenericSparsityPattern & sparsity_pattern ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::copy() const -> Matrix *
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::size( size_t dim ) const -> size_t
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::zero() -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::apply( FinalizeType final ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::disp( size_t precision ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t M, size_t N ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t M, size_t N, bool distributed ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::get( real *         block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::set( const real *   block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::add( const real *   block,
                  size_t         m,
                  const size_t * rows,
                  size_t         n,
                  const size_t * cols ) -> void
{
  // FIXME
}
//-----------------------------------------------------------------------------

auto Matrix::norm( std::string norm_type ) const -> real
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::getrow( size_t                  row,
                     std::vector< size_t > & columns,
                     std::vector< real > &   values ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::setrow( size_t                        row,
                     const std::vector< size_t > & columns,
                     const std::vector< real > &   values ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::zero( size_t m, const size_t * rows ) -> void
{
  // FIXME
}
//-----------------------------------------------------------------------------

auto Matrix::ident( size_t m, const size_t * rows ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::dup( GenericMatrix & A ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::mult( const GenericVector & x,
                   GenericVector &       y,
                   bool                  transposed ) const -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::operator*=( real a ) -> const Matrix &
{
  // FIXME
}
//-----------------------------------------------------------------------------

auto Matrix::operator/=( real a ) -> const Matrix &
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::operator=( const GenericMatrix & A ) -> const GenericMatrix &
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::nz() const -> size_t
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::factory() const -> LinearAlgebraFactory &
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::norm( const Norm type ) const -> real
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::operator=( const Matrix & A ) -> const Matrix &
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::operator+=( const Matrix & A ) -> const Matrix &
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t M, size_t N, std::vector< size_t > const & nz )
  -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::init( size_t                        M,
                   size_t                        N,
                   std::vector< size_t > const & d_nzrow,
                   std::vector< size_t > const & o_nzrow ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

auto Matrix::getrows_offproc( _ordered_set< size_t > const & rows ) -> void
{
  // FIXME
}

//-----------------------------------------------------------------------------

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS
