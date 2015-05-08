// Copyright (C) 2007 Ola Skavhaug
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrom 2008.
//
// First added:  2007-11-30
// Last changed: 2008-01-24

#ifndef __GENERIC_SPARSITY_PATTERN_H
#define __GENERIC_SPARSITY_PATTERN_H

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
    GenericSparsityPattern() {}

    /// Destructor
    virtual ~GenericSparsityPattern() {};

    /// Initialise sparsity pattern for a matrix with total number of rows and
    /// columns and set if it is distributed
    virtual void init(uint rank, uint const * dims, bool distributed) = 0;

    /// Initialize sparsity pattern for a generic tensor
    virtual void init(uint rank, uint const * dims) = 0;

    /// Initialise sparsity pattern for a parallel generic tensor
    /// [TODO: Deprecate]
    virtual void pinit(uint rank, uint const * dims) = 0;

    /// Insert non-zero entry
    virtual void insert(uint const * num_rows, uint const * const * rows) = 0;

    /// Insert non-zero entry [TODO: Deprecate]
    virtual void pinsert(uint const * num_rows, uint const * const * rows) = 0;

    /// Return local size
    virtual uint size(uint n) const = 0;

    /// Finalize sparsity pattern (needed by most parallel la backends)
    virtual void apply() = 0;

    /// Clear
    virtual void clear() = 0;

    /// Is blocked
    virtual bool is_blocked() const = 0;

    /// Is distributed
    virtual bool is_distributed() const = 0;

    /// Initialize process range
    virtual void initRange(uint num_local) = 0;

    /// Return array with number of non-zeroes per row
    virtual void numNonZeroPerRow(uint nzrow[]) const = 0;

    /// Return array with number of non-zeroes per row diagonal and offdiagonal
    /// for process_number
    virtual void numNonZeroPerRow(uint process_number, uint d_nzrow[],
                                  uint o_nzrow[]) const = 0;

    /// Return total number of non-zeroes
    virtual uint numNonZero() const = 0;

    /// Return the pattern
    virtual Array< _set<int> > const& pattern() const = 0;

    /// Display sparsity pattern
    virtual void disp() const = 0;

  };

}

#endif
