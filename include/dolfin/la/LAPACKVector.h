// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LAPACK_VECTOR_H
#define __DOLFIN_LAPACK_VECTOR_H

#include <dolfin/common/Variable.h>
#include <dolfin/common/types.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// This class provides a simple wrapper for matrix data for use
/// with LAPACK (column-major ordering).
///
/// This class does currently not implement the GenericVector
/// interface but may possibly be extended to do so in the future.

class LAPACKVector : public Variable
{
public:
  /// Create M x N matrix
  LAPACKVector( size_t M );

  /// Destructor
  ~LAPACKVector();

  /// Return size of vector
  auto size() const -> size_t
  {
    return M_;
  }

  auto values() -> double *
  {
    return values_;
  }

  /// Access entry i
  auto operator[]( size_t i ) -> double &
  {
    return values_[i];
  }

  /// Access entry i, const version
  auto operator[]( size_t i ) const -> double
  {
    return values_[i];
  }

  /// Return informal string representation (pretty-print)
  auto str( bool verbose ) const -> std::string;

private:
  // Number of rows and columns
  size_t M_;

  // Values, stored column-major
  double * values_;
};

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
