// Copyright (C) 2007 Garth N. Wells
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SPARSITY_PATTERN_H
#define __DOLFIN_SPARSITY_PATTERN_H

#include <dolfin/la/GenericSparsityPattern.h>

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

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

  /// Create sparsity pattern for given global dimensions and local ranges.
  /// If range is a nullptr pointer the pattern is assumed to be serial.
  SparsityPattern(uint rank, uint const * dim, uint const * range = nullptr);

  /// Destructor
  ~SparsityPattern();

  //--- INTERFACE -------------------------------------------------------------

  /// Initialize with given tensor rank, global dimensions and local ranges.
  /// If range is a nullptr pointer the pattern is assumed to be serial
  void init(uint rank, uint const * dim, uint const * range = nullptr);

  /// Clear
  void clear();

  /// Insert non-zero entries
  void insert(uint const * num, uint const * const * idx);

  /// Return local size for given dimension
  uint size(uint i) const;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  void apply();

  /// Is blocked
  bool is_blocked() const;

  /// Is distributed
  bool is_distributed() const;

  /// Return array with number of non-zeroes per local row
  void numNonZeroPerRow(uint nzrow[]) const;

  /// Return array with number of non-zeroes per row for the given process rank
  /// and split between entries in the diagonal and off-diagonal portion of the
  /// matrix
  void numNonZeroPerRow(uint p_rank, uint d_nzrow[], uint o_nzrow[]) const;

  /// Return total number of non-zeroes
  uint numNonZero() const;

  /// Display sparsity pattern
  void disp() const;

  //---------------------------------------------------------------------------

  /// Return array with row range for process rank
  void get_range(uint p_rank, uint range[]);

  /// Return number of local rows for process rank
  uint range_size(uint p_rank) const;

  ///
  void set_blocked();

private:

  /// Tensor rank
  uint rank_;

  /// Dimensions
  uint * dim_;

  /// Range -array of size + 1 where size is size + 1:
  ///    range[rank], range[rank+1] is the range for processor
  uint ** range_;

  /// Direct access to local range
  uint ** local_range_;

  /// Flags
  bool initialized_;
  bool finalized_;
  bool blocked_;
  bool distributed_;

  /// Sparsity pattern represented as an array of sets.
  /// Each set corresponds to a row in the local range and contains the column
  /// positions of nonzero entries.

  /// Diagonal portion: submatrix such that row and column
  /// indices are in-range
  _ordered_set<uint> * d_entries_;
  uint d_count_;

  /// Off-diagonal portion: entries such that only column indices are off-range
  _ordered_set<uint> * o_entries_;
  uint o_count_;

  /// Additionally provide data structure to store remote entries i,e such that
  /// row indices are not in-range
  _ordered_map<uint, _ordered_set<uint> > r_entries_;

};

//-----------------------------------------------------------------------------
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
inline uint SparsityPattern::size( uint i ) const
{
  return ( local_range_[i][1] - local_range_[i][0] );
}
//-----------------------------------------------------------------------------
inline void SparsityPattern::get_range( uint p_rank, uint range[] )
{
  dolfin_assert( distributed_ );
  // For a serial pattern p_rank is only zero
  std::copy( &range_[0][p_rank], &range_[0][p_rank + 1], range );
}
//-----------------------------------------------------------------------------
inline uint SparsityPattern::range_size( uint p_rank ) const
{
  dolfin_assert( distributed_ );
  return range_[0][p_rank + 1] - range_[0][p_rank];
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_SPARSITY_PATTERN_H */
