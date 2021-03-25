// Copyright (C) 2020 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_MATRIX_H
#define __DOLFIN_TRILINOS_MATRIX_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/trilinos/TrilinosObject.h>

#include <Tpetra_CrsMatrix.hpp>

namespace dolfin
{

namespace trilinos
{

class Matrix : public GenericMatrix, public Object, public Variable
{
public:
  using LO = int;
  using GO = size_t;

  static constexpr GO const indexBase = 0;

	  /// Matrix type (scalar, local index, global index)
  using TPMatrix = Tpetra::CrsMatrix< real, LO, GO >;

  /// Graph type (local index, global index)
  using TPGraph  = Tpetra::CrsGraph< LO, GO >;

  /// Map type (local index, global index)
  using TPMap    = Tpetra::Map< LO, GO >;

public:
  enum Norm
  {
    l1,
    linf,
    frobenius
  };

public:
  /// Create empty matrix
  explicit Matrix();

  /// Create matrix of local dimension M x N
  Matrix( size_t M, size_t N, bool distributed = true );

  /// Copy constructor
  explicit Matrix( const Matrix & A );

  // /// Create matrix from given Trilinos Mat pointer
  // explicit Matrix( Mat A );

  /// Destructor
  ~Matrix() override;

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  auto copy() const -> Matrix * override;

  /// Return size of given dimension
  auto size( size_t dim ) const -> size_t override;

  /// Set all entries to zero and keep any sparse structure
  auto zero() -> void override;

  /// Finalize assembly of tensor
  auto apply( FinalizeType finaltype = FINALIZE ) -> void override;

  /// Display tensor
  auto disp( size_t precision = 0 ) const -> void override;

  /// Initialize zero tensor using sparsity pattern
  auto init( const GenericSparsityPattern & sparsity_pattern ) -> void override;

  //--- Implementation of the GenericMatrix interface --

  /// Initialize matrix of local dimension M x N, distributed by default
  auto init( size_t M, size_t N ) -> void override;

  /// Initialize matrix of local dimension M x N, distributed if specified
  auto init( size_t M, size_t N, bool distributed ) -> void override;

  /// Get block of values
  auto get( real *         block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) const -> void override;

  /// Set block of values
  auto set( const real *   block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) -> void override;

  /// Add block of values
  auto add( const real *   block,
            size_t         m,
            const size_t * rows,
            size_t         n,
            const size_t * cols ) -> void override;

  /// Return norm of matrix
  auto norm( std::string norm_type = "frobenius" ) const -> real override;

  /// Get non-zero values of given row
  auto getrow( size_t                  row,
               std::vector< size_t > & columns,
               std::vector< real > &   values ) const -> void override;

  /// Set values for given row
  auto setrow( size_t                        row,
               const std::vector< size_t > & columns,
               const std::vector< real > &   values ) -> void override;

  /// Set given rows to zero
  auto zero( size_t m, const size_t * rows ) -> void override;

  /// Set given rows to identity matrix
  auto ident( size_t m, const size_t * rows ) -> void override;

  // Matrix-vector product, y = Ax
  auto mult( const GenericVector & x,
             GenericVector &       y,
             bool                  transposed = false ) const -> void override;

  /// Multiply matrix by given number
  auto operator*=( real a ) -> const Matrix & override;

  /// Divide matrix by given number
  auto operator/=( real a ) -> const Matrix & override;

  /// Assignment operator
  auto operator=( const GenericMatrix & A ) -> const GenericMatrix & override;

  /// Get number of non-zeros in the matrix
  auto nz() const -> size_t override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  auto factory() const -> LinearAlgebraFactory & override;

  //--- Special Functions ---

  /// Return Trilinos Mat pointer
  // auto mat() const -> Mat;

  /// Return norm of matrix
  auto norm( const Norm type = l1 ) const -> real;

  /// Assignment operator
  auto operator=( const Matrix & A ) -> const Matrix &;

  /// Matrix axpy, Y = a X+ Y
  auto operator+=( const Matrix & A ) -> const Matrix &;

private:
  //
  auto clear() -> void;

  // Print info
  // auto print( MatInfo const & info ) const -> void;

  // The matrix
  Teuchos::RCP< TPMatrix > mat_;

  // // Sub-matrices
  // Mat * AA_sub { nullptr };

  // // True if the matrix is distributed
  // bool is_distributed_ { false };

  // PetscInt rstart_ { 0 };
  // PetscInt rend_ { 0 };

  // _map< int, int > mapping_;
};

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_MATRIX_H
