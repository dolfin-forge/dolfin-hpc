// Copyright (C) 2007 Ola Skavhaug
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrom 2008.
//
// First added:  2007-11-30
// Last changed: 2008-01-24

#ifndef __DOLFIN_GENERIC_SPARSITY_PATTERN_H
#define __DOLFIN_GENERIC_SPARSITY_PATTERN_H

#include <dolfin/common/types.h>

namespace dolfin
{

template<class T> class Array;

/// Base class for sparsity patterns of vectors/matrices. Concrete
/// sub classes can be used to initialize vectors and sparse
/// matrices.

class GenericSparsityPattern
{

public:

  /// Constructor
  GenericSparsityPattern()
  {
  }

  /// Destructor
  virtual ~GenericSparsityPattern()
  {
  }

  /// Initialize sparsity pattern for a generic tensor
  virtual void init(uint rank, uint const * dim, uint const * range) = 0;

  /// Insert non-zero entry
  virtual void insert(uint const * num_rows, uint const * const * rows) = 0;

  /// Return local size
  virtual uint size(uint i) const = 0;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  virtual void apply() = 0;

  /// Clear
  virtual void clear() = 0;

  /// Is blocked
  virtual bool is_blocked() const = 0;

  /// Is distributed
  virtual bool is_distributed() const = 0;

  /// Return array with number of non-zeroes per row
  virtual void numNonZeroPerRow(uint nzrow[]) const = 0;

  /// Return array with number of non-zeroes per row split between
  /// diagonal and off-diagonal for process rank
  virtual void numNonZeroPerRow(uint p_rank, uint d_nzrow[],
                                uint o_nzrow[]) const = 0;

  /// Return total number of non-zeroes
  virtual uint numNonZero() const = 0;

  /// Display sparsity pattern
  virtual void disp() const = 0;

};

} /* namespace dolfin */

#endif /* __DOLFIN_GENERIC_SPARSITY_PATTERN_H */
