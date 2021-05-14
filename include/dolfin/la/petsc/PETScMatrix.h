// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_MATRIX_H
#define __DOLFIN_PETSC_MATRIX_H

#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericMatrix.h>

#include <vector>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScObject.h>

#include <petscmat.h>

namespace dolfin
{

class PETScVector;

/// This class provides a simple matrix class based on PETSc.
/// It is a simple wrapper for a PETSc matrix pointer (Mat)
/// implementing the GenericMatrix interface.
///
/// The interface is intentionally simple. For advanced usage,
/// access the PETSc Mat pointer using the function mat() and
/// use the standard PETSc interface.

class PETScMatrix : public GenericMatrix, public PETScObject, public Variable
{

public:
  enum Norm
  {
    l1,
    linf,
    frobenius
  };

  /// Create empty matrix
  explicit PETScMatrix();

  /// Create matrix of local dimension M x N
  PETScMatrix( size_t M, size_t N, bool distributed = true );

  /// Copy constructor
  explicit PETScMatrix( const PETScMatrix & A );

  /// Create matrix from given PETSc Mat pointer
  explicit PETScMatrix( Mat A );

  /// Destructor
  ~PETScMatrix() override;

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  void init( const GenericSparsityPattern & sparsity_pattern ) override;

  /// Return copy of tensor
  auto copy() const -> PETScMatrix * override;

  /// Return size of given dimension
  auto size( size_t dim ) const -> size_t override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType final = FINALIZE ) override;

  /// Display tensor
  void disp( size_t precision = 0 ) const override;

  //--- Implementation of the GenericMatrix interface --

  /// Initialize matrix of local dimension M x N, distributed by default
  void init( size_t M, size_t N ) override;

  /// Initialize matrix of local dimension M x N, distributed if specified
  void init( size_t M, size_t N, bool distributed ) override;

  /// Get block of values
  void get( real *         block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) const override;

  /// Set block of values
  void set( const real *   block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) override;

  /// Add block of values
  void add( const real *   block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) override;

  /// Return norm of matrix
  auto norm( std::string norm_type = "frobenius" ) const -> real override;

  /// Get non-zero values of given row
  void getrow( size_t                  row,
               std::vector< size_t > & columns,
               std::vector< real > &   values ) const override;

  /// Set values for given row
  void setrow( size_t                        row,
               const std::vector< size_t > & columns,
               const std::vector< real > &   values ) override;

  /// Set given rows to zero
  void zero( size_t m, const size_t * rows ) override;

  /// Set given rows to identity matrix
  void ident( size_t m, const size_t * rows ) override;

  /// Duplicate matrix
  void dup( GenericMatrix & A );

  // Matrix-vector product, y = Ax
  void mult( const GenericVector & x,
             GenericVector &       y,
             bool                  transposed = false ) const override;

  /// Multiply matrix by given number
  auto operator*=( real a ) -> const PETScMatrix & override;

  /// Divide matrix by given number
  auto operator/=( real a ) -> const PETScMatrix & override;

  /// Assignment operator
  auto operator=( const GenericMatrix & A ) -> const GenericMatrix & override;

  /// Get number of non-zeros in the matrix
  auto nz() const -> size_t override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  auto factory() const -> LinearAlgebraFactory & override;

  //--- Special PETScFunctions ---

  /// Return PETSc Mat pointer
  auto mat() const -> Mat;

  /// Return norm of matrix
  auto norm( const Norm type = l1 ) const -> real;

  /// Assignment operator
  auto operator=( const PETScMatrix & A ) -> const PETScMatrix &;

  /// Matrix axpy, Y = a X+ Y
  auto operator+=( const PETScMatrix & A ) -> const PETScMatrix &;

private:
  //
  void clear();

  // Initialize M x N matrix with a given number of nonzeros per row
  void init( size_t M, size_t N, std::vector< size_t > const & nz );

  // Initialize M x N matrix with a given number of nonzeros per row diagonal
  // and off-diagonal
  void init( size_t M, size_t N,
             std::vector< size_t > const & d_nzrow,
             std::vector< size_t > const & o_nzrow );

  ///
  void getrows_offproc( _ordered_set< size_t > const & rows );

  // Print info
  void print( MatInfo const & info ) const;

  // Matrix
  Mat A { nullptr };

  // Sub-matrices
  Mat * AA_sub { nullptr };

  // True if the matrix is distributed
  bool is_distributed_ { false };

  PetscInt rstart_ { 0 };
  PetscInt rend_ { 0 };

  _map< int, int > mapping_;
};

//-----------------------------------------------------------------------------
inline void PETScMatrix::init( size_t M, size_t N )
{
  init( M, N, true );
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::size( size_t dim ) const -> size_t
{
  PetscInt M = 0;
  PetscInt N = 0;
  MatGetSize( A, &M, &N );
  return ( dim == 0 ? M : N );
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::nz() const -> size_t
{
  MatInfo info;
  MatGetInfo( A, MAT_GLOBAL_SUM, &info );

  return info.nz_used;
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::get( real *         block,
                              size_t         m,
                              size_t const * rows,
                              size_t         n,
                              size_t const * cols ) const
{
  dolfin_assert( A );

  PetscInt M_ = m;
  PetscInt N_ = n;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );
  std::vector< PetscInt > cols_( cols, cols + n );

  MatGetValues( A, M_, rows_.data(), N_, cols_.data(), block );
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::set( real const *   block,
                              size_t         m,
                              size_t const * rows,
                              size_t         n,
                              size_t const * cols )
{
  dolfin_assert( A );

  PetscInt M_ = m;
  PetscInt N_ = n;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );
  std::vector< PetscInt > cols_( cols, cols + n );

  MatSetValues( A, M_, rows_.data(), N_, cols_.data(), block, INSERT_VALUES );
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::add( real const *   block,
                              size_t         m,
                              size_t const * rows,
                              size_t         n,
                              size_t const * cols )
{
  dolfin_assert( A );

  PetscInt M_ = m;
  PetscInt N_ = n;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );
  std::vector< PetscInt > cols_( cols, cols + n );

  MatSetValues( A, M_, rows_.data(), N_, cols_.data(), block, ADD_VALUES );
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::setrow( size_t                        row,
                                 const std::vector< size_t > & columns,
                                 const std::vector< real > &   values )
{
  set( values.data(), 1, &row, columns.size(), columns.data() );
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::apply( FinalizeType finaltype )
{
  if ( finaltype == FINALIZE )
  {
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
  }
  else if ( finaltype == FLUSH )
  {
    MatAssemblyBegin( A, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd( A, MAT_FLUSH_ASSEMBLY );
  }
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::zero()
{
  MatZeroEntries( A );
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::operator+=( const PETScMatrix & A )
  -> const PETScMatrix &
{
  dolfin_assert( this->A );
  MatAXPY( this->A, 1.0, A.A, SAME_NONZERO_PATTERN );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::operator*=( real a ) -> const PETScMatrix &
{
  dolfin_assert( A );
  MatScale( A, a );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::operator/=( real a ) -> const PETScMatrix &
{
  dolfin_assert( A );
  MatScale( A, 1.0 / a );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::operator=( const GenericMatrix & A )
  -> const GenericMatrix &
{
  if ( &A != this )
  {
    MatCopy( A.down_cast< PETScMatrix >().A, this->A, SAME_NONZERO_PATTERN );
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::operator=( const PETScMatrix & A )
  -> const PETScMatrix &
{
  if ( &A != this )
  {
    MatCopy( A.A, ( this->A ), SAME_NONZERO_PATTERN );
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::dup( GenericMatrix & A )
{
  MatDuplicate( A.down_cast< PETScMatrix >().A, MAT_COPY_VALUES, &this->A );
}
//-----------------------------------------------------------------------------
inline auto PETScMatrix::mat() const -> Mat
{
  return A;
}

} /* namespace dolfin */

#endif /* HAVE_PETSC */

#endif /* __DOLFIN_PETSC_MATRIX_H */
