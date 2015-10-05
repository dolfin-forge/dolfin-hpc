// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-12-14
// Last changed: 2009-12-15

#ifndef __DOLFIN_LAPACK_VECTOR_H
#define __DOLFIN_LAPACK_VECTOR_H

#include <dolfin/common/types.h>
#include <dolfin/common/Variable.h>

#include <cstring>
#include <string>

namespace dolfin
{

/// This class provides a simple wrapper for matrix data for use
/// with LAPACK (column-major ordering).
///
/// This class does currently not implement the GenericVector
/// interface but may possibly be extended to do so in the future.

class LAPACKVector : public Variable
{
public:

  /// Create M x N matrix
  LAPACKVector(uint M);

  /// Destructor
  ~LAPACKVector();

  /// Return size of vector
  uint size() const
  {
    return M_;
  }

  double * values()
  {
    return values_;
  }

  /// Access entry i
  double& operator[](uint i)
  {
    return values_[i];
  }

  /// Access entry i, const version
  double operator[](uint i) const
  {
    return values_[i];
  }

  /// Return informal string representation (pretty-print)
  std::string str(bool verbose) const;

private:

  // Number of rows and columns
  uint M_;

  // Values, stored column-major
  double * values_;

};

}

#endif
