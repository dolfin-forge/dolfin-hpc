// Copyright (C) 2007 Garth N. Wells
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2007-2008.
// Modified by Aurélien Larcher, 2015.
//
// First added:  2007-03-13
// Last changed: 2008-05-15

#ifndef __SPARSITY_PATTERN_H
#define __SPARSITY_PATTERN_H

#include "GenericSparsityPattern.h"

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

#include <set>

namespace dolfin
{

/// This class represents the sparsity pattern of a vector/matrix. It can be
/// used to initalise vectors and sparse matrices. It must be initialised
/// before use.

class SparsityPattern: public GenericSparsityPattern
{
public:

  /// Create empty sparsity pattern
  SparsityPattern();

  /// Create sparsity pattern for tensor of given dimensions
  SparsityPattern(uint rank, uint const * dims, bool distributed);

  /// Destructor
  ~SparsityPattern();

  //--- INTERFACE -------------------------------------------------------------

  /// Initialise sparsity pattern for a matrix with total number of rows and
  /// columns and set if it is distributed
  void init(uint rank, uint const * dims, bool distributed);

  /// Initialise sparsity pattern for a matrix with total number of rows and
  /// columns
  void init(uint rank, uint const * dims);

  /// Initialise sparsity pattern for a parallel matrix with total number of
  /// rows and columns
  void pinit(uint rank, uint const * dims);

  /// Insert non-zero entries
  void insert(uint const * num_rows, uint const * const * rows);

  /// Insert non-zero entry for parallel matrices
  void pinsert(uint const * num_rows, uint const * const * rows);

  /// Return local size
  uint size(uint n) const;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  void apply();

  /// Return underlying sparsity pattern after apply has been called
  Array< _set<int> > const& pattern() const;

  /// Clear
  void clear();

  /// Is blocked
  bool is_blocked() const;

  /// Is distributed
  bool is_distributed() const;

  /// Initialize process range
 void initRange(uint num_local);

  /// Return array with number of non-zeroes per row
  void numNonZeroPerRow(uint nzrow[]) const;

  /// Return array with number of non-zeroes per row diagonal and offdiagonal
  /// for process_number
  void numNonZeroPerRow(uint process_number, uint d_nzrow[],
                        uint o_nzrow[]) const;

  /// Return total number of non-zeroes
  uint numNonZero() const;

  /// Display sparsity pattern
  void disp() const;

  //---------------------------------------------------------------------------

  /// Return array with row range for process_number
  void processRange(uint process_number, uint local_range[]);

  /// Return number of local rows for process_number
  uint numLocalRows(uint process_number) const;

  ///
  void set_blocked();

private:

  /// Sparsity pattern represented as an vector of sets. Each set corresponds
  /// to a row, and the set contains the column positions of nonzero entries
  /// When run in parallel this vector contains diagonal non-zeroes
#if __SUNPRO_CC
  std::map<uint,  std::set<int> > sparsity_pattern;
#else
  std::map<uint const,  std::set<int> > sparsity_pattern;
#endif

  /// Sparsity pattern for off diagonal represented as vector of sets. Each
  /// set corresponds to a row, and the set contains the column positions of nonzero entries
#if __SUNPRO_CC
  std::map<uint,  std::set<int> > o_sparsity_pattern;
#else
  std::map<uint const,  std::set<int> > o_sparsity_pattern;
#endif

  // Dimensions
  uint dim[2];

  //range -array of size + 1 where size is numProcesses + 1.
  //range[rank], range[rank+1] is the range for processor
  uint * range;

  std::vector<int> off_processor;

  bool initialized_;
  bool finalized_;
  bool blocked_;
  bool distributed_;

  // Dummy return value
  Array< _set<int> > pattern_;
};

//--- INLINES -----------------------------------------------------------------

inline void SparsityPattern::set_blocked()
{
  blocked_ = true;
}

//-----------------------------------------------------------------------------
inline bool SparsityPattern::is_blocked() const
{
  return blocked_;
}

//-----------------------------------------------------------------------------
inline bool SparsityPattern::is_distributed() const
{
  return distributed_;
}

//-----------------------------------------------------------------------------

}
#endif
