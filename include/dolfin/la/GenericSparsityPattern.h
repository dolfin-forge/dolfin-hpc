// Copyright (C) 2007 Ola Skavhaug
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_SPARSITY_PATTERN_H
#define __DOLFIN_GENERIC_SPARSITY_PATTERN_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// Base class for sparsity patterns of vectors/matrices. Concrete
/// sub classes can be used to initialize vectors and sparse
/// matrices.

class GenericSparsityPattern
{

public:

  /// Constructor
  GenericSparsityPattern() = default;

  /// Destructor
  virtual ~GenericSparsityPattern() = default;

  /// Initialize sparsity pattern for a generic tensor
  virtual void init(size_t rank, size_t const * dim, size_t const * range) = 0;

  /// Insert non-zero entry
  virtual void insert(size_t const * num_rows, size_t const * const * rows) = 0;

  /// Return local size
  virtual auto size(size_t i) const -> size_t = 0;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  virtual void apply() = 0;

  /// Clear
  virtual void clear() = 0;

  /// Is blocked
  virtual auto is_blocked() const -> bool = 0;

  /// Is distributed
  virtual auto is_distributed() const -> bool = 0;

  /// Return array with number of non-zeroes per row
  virtual void numNonZeroPerRow(size_t nzrow[]) const = 0;

  /// Return array with number of non-zeroes per row split between
  /// diagonal and off-diagonal for process rank
  virtual void numNonZeroPerRow(size_t p_rank, size_t d_nzrow[],
                                size_t o_nzrow[]) const = 0;

  /// Return total number of non-zeroes
  virtual auto numNonZero() const -> size_t = 0;

  /// Display sparsity pattern
  virtual void disp() const = 0;

};

} /* namespace dolfin */

#endif /* __DOLFIN_GENERIC_SPARSITY_PATTERN_H */
