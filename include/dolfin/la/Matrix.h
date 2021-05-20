// Copyright (C) 2006-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MATRIX_H
#define __DOLFIN_MATRIX_H

#include <dolfin/common/Variable.h>
#include <dolfin/la/DefaultFactory.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/main/PE.h>

#include <fstream>
#include <sstream>

namespace dolfin
{

/// This class provides the default DOLFIN matrix class,
/// based on the default DOLFIN linear algebra backend.

class Matrix : public GenericMatrix, public Variable
{
public:
  /// Create empty matrix
  Matrix()
    : Variable( "A", "DOLFIN matrix" )
  {
    DefaultFactory factory;
    matrix_ = factory.createMatrix();
  }

  /// Create M x N matrix distributed by default
  Matrix( size_t M, size_t N )
    : Variable( "A", "DOLFIN matrix" )
    , matrix_( nullptr )
  {
    DefaultFactory factory;
    matrix_ = factory.createMatrix();
    matrix_->init( M, N );
  }

  /// Create M x N matrix distributed if specified
  Matrix( size_t M, size_t N, bool distributed )
    : Variable( "A", "DOLFIN matrix" )
    , matrix_( nullptr )
  {
    DefaultFactory factory;
    matrix_ = factory.createMatrix();
    matrix_->init( M, N, distributed );
  }

  /// Copy constructor
  explicit Matrix( const Matrix & A )
    : Variable( "A", "DOLFIN matrix" )
    , matrix_( A.matrix_->copy() )
  {
  }

  /// Destructor
  ~Matrix() override
  {
    delete matrix_;
  }

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  /// The tensor is distributed if the sparsity pattern is distributed
  void init( const GenericSparsityPattern & sparsity_pattern ) override;

  /// Return copy of tensor
  Matrix * copy() const override;

  /// Return size of given dimension
  size_t size( size_t dim ) const override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE ) override;

  /// Display tensor
  void disp( size_t precision = 2 ) const override;

  //--- Implementation of the GenericMatrix interface ---

  /// Initialize M x N matrix and distribute by default
  void init( size_t M, size_t N ) override;

  /// Initialize M x N matrix and distribute if specified
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
  double norm( std::string norm_type = "frobenius" ) const override;

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

  // Matrix-vector product, y = Ax
  void mult( const GenericVector & x,
             GenericVector &       y,
             bool                  transposed = false ) const override;

  /// Multiply matrix by given number
  const Matrix & operator*=( real a ) override;

  /// Divide matrix by given number
  const Matrix & operator/=( real a ) override;

  /// Assignment operator
  const GenericMatrix & operator=( const GenericMatrix & A ) override;

  /// Get number of non-zeros in the matrix
  size_t nz() const override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory & factory() const override;

  //--- Special functions, intended for library use only ---

  /// Return concrete instance / unwrap (const)
  const GenericMatrix * instance() const override;

  /// Return concrete instance / unwrap (non-const version)
  GenericMatrix * instance() override;

  //--- Special Matrix functions ---

  /// Assignment operator
  const Matrix & operator=( const Matrix & A );

  void spy() const;

private:
  // Pointer to concrete implementation
  GenericMatrix * matrix_ { nullptr };
};

//-----------------------------------------------------------------------------

inline void Matrix::init( const GenericSparsityPattern & sparsity_pattern )
{
  matrix_->init( sparsity_pattern );
}

//-----------------------------------------------------------------------------

inline Matrix * Matrix::copy() const
{
  Matrix * A = new Matrix();
  delete A->matrix_;
  A->matrix_ = matrix_->copy();
  return A;
}

//-----------------------------------------------------------------------------

inline size_t Matrix::size( size_t dim ) const
{
  return matrix_->size( dim );
}

//-----------------------------------------------------------------------------

inline void Matrix::zero()
{
  matrix_->zero();
}

//-----------------------------------------------------------------------------

inline void Matrix::apply( FinalizeType finaltype )
{
  matrix_->apply( finaltype );
}

//-----------------------------------------------------------------------------

inline void Matrix::disp( size_t precision ) const
{
  matrix_->disp( precision );
}

//-----------------------------------------------------------------------------

inline void Matrix::init( size_t M, size_t N )
{
  matrix_->init( M, N );
}

//-----------------------------------------------------------------------------

inline void Matrix::init( size_t M, size_t N, bool distributed )
{
  matrix_->init( M, N, distributed );
}

//-----------------------------------------------------------------------------

inline void Matrix::get( real *         block,
                         size_t         m,
                         const size_t * rows,
                         size_t         n,
                         const size_t * cols ) const
{
  matrix_->get( block, m, rows, n, cols );
}

//-----------------------------------------------------------------------------

inline void Matrix::set( const real *   block,
                         size_t         m,
                         const size_t * rows,
                         size_t         n,
                         const size_t * cols )
{
  matrix_->set( block, m, rows, n, cols );
}

//-----------------------------------------------------------------------------

inline void Matrix::add( const real *   block,
                         size_t         m,
                         const size_t * rows,
                         size_t         n,
                         const size_t * cols )
{
  matrix_->add( block, m, rows, n, cols );
}

//-----------------------------------------------------------------------------

inline double Matrix::norm( std::string norm_type ) const
{
  return matrix_->norm( norm_type );
}

//-----------------------------------------------------------------------------

inline void Matrix::getrow( size_t                  row,
                            std::vector< size_t > & columns,
                            std::vector< real > &   values ) const
{
  matrix_->getrow( row, columns, values );
}

//-----------------------------------------------------------------------------

inline void Matrix::setrow( size_t                        row,
                            const std::vector< size_t > & columns,
                            const std::vector< real > &   values )
{
  matrix_->setrow( row, columns, values );
}

//-----------------------------------------------------------------------------

inline void Matrix::zero( size_t m, const size_t * rows )
{
  matrix_->zero( m, rows );
}

//-----------------------------------------------------------------------------

inline void Matrix::ident( size_t m, const size_t * rows )
{
  matrix_->ident( m, rows );
}

//-----------------------------------------------------------------------------

inline void Matrix::mult( const GenericVector & x,
                          GenericVector &       y,
                          bool                  transposed ) const
{
  matrix_->mult( x, y, transposed );
}

//-----------------------------------------------------------------------------

inline const Matrix & Matrix::operator*=( real a )
{
  *matrix_ *= a;
  return *this;
}
//-----------------------------------------------------------------------------

inline const Matrix & Matrix::operator/=( real a )
{
  *matrix_ /= a;
  return *this;
}
inline const GenericMatrix & Matrix::operator=( const GenericMatrix & A )
{
  *matrix_ = A;
  return *this;
}

//-----------------------------------------------------------------------------

inline size_t Matrix::nz() const
{
  return matrix_->nz();
}

//-----------------------------------------------------------------------------

inline LinearAlgebraFactory & Matrix::factory() const
{
  return matrix_->factory();
}

//-----------------------------------------------------------------------------

inline const GenericMatrix * Matrix::instance() const
{
  return matrix_;
}

//-----------------------------------------------------------------------------

inline GenericMatrix * Matrix::instance()
{
  return matrix_;
}

//-----------------------------------------------------------------------------

inline const Matrix & Matrix::operator=( const Matrix & A )
{
  *matrix_ = *A.matrix_;
  return *this;
}

//-----------------------------------------------------------------------------

inline void Matrix::spy() const
{
  if ( PE::size() == 1 )
  {
    std::stringstream ss;
    ss << "A" << matrix_->size( 0 ) * matrix_->size( 1 ) << ".xpm" << std::ends;
    std::ofstream     Afile( ss.str().c_str() );
    real              val = 0.;
    std::stringstream xpm;
    xpm << "/* XPM */\n"
        << "static char * map_xpm = {\n"
        << "/* width height number_of_colors chars_per_pixel */\n"
        << "\"" << matrix_->size( 0 ) << " " << matrix_->size( 1 ) << " 2 1\",\n"
        << "/* intensity levels */\n"
        << "\"0 c #ffffff\",\n"
        << "\"1 c none\",\n"
        << "/* map */\n";
    Afile << xpm.str();
    for ( size_t i = 0; i < matrix_->size( 0 ); ++i )
    {
      std::stringstream row;
      row << "\"";
      for ( size_t j = 0; j < matrix_->size( 1 ); ++j )
      {
        val = matrix_->getitem( std::pair< size_t, size_t >( i, j ) );
        if ( val > 1.0e-14 )
        {
          row << "1";
        }
        else
        {
          row << "0";
        }
      }
      row << "\",\n";
      Afile << row.str();
    }
    Afile << "};";
  }
}

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
