// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LAPACK_MATRIX_H
#define __DOLFIN_LAPACK_MATRIX_H

#include <dolfin/common/types.h>
#include <dolfin/common/Variable.h>
#include <dolfin/log/log.h>

#include <cstring>
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
  uint size(uint dim) const
  {
    dolfin_assert(dim < 2);
    return (dim == 0 ? M_ : N_);
  }

  double * values()
  {
    return values_;
  }

  /// Access entry (i, j)
  double& operator()(uint i, uint j)
  {
    return values_[j * M_ + i];
  }

  /// Access entry (i, j), const version
  double operator()(uint i, uint j) const
  {
    return values_[j * M_ + i];
  }

  /// Return informal string representation (pretty-print)
  std::string str(bool verbose) const;

private:

  // Number of rows and columns
  uint M_;
  uint N_;

  // Values, stored column-major
  double * values_;

};

}

#endif
