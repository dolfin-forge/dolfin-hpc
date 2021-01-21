// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LAPACK_MATRIX_H
#define __DOLFIN_LAPACK_MATRIX_H

#include <dolfin/common/assert.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Variable.h>

#include <string>

namespace dolfin
{

/// This class provides a simple wrapper for matrix data for use
/// with LAPACK (column-major ordering).
///
/// This class does currently not implement the GenericMatrix
/// interface but may possibly be extended to do so in the future.

class LAPACKMatrix : public Variable
{
public:

  /// Create M x N matrix
  LAPACKMatrix(uint M, uint N);

  /// Destructor
  ~LAPACKMatrix();

  /// Return size of given dimension
  auto size(uint dim) const -> uint
  {
    dolfin_assert(dim < 2);
    return (dim == 0 ? M_ : N_);
  }

  auto values() -> double *
  {
    return values_;
  }

  /// Access entry (i, j)
  auto operator()(uint i, uint j) -> double&
  {
    return values_[j * M_ + i];
  }

  /// Access entry (i, j), const version
  auto operator()(uint i, uint j) const -> double
  {
    return values_[j * M_ + i];
  }

  /// Return informal string representation (pretty-print)
  auto str(bool verbose) const -> std::string;

private:

  // Number of rows and columns
  uint M_;
  uint N_;

  // Values, stored column-major
  double * values_;

};

}

#endif
