// Copyright (C) 2007 Garth N. Wells
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SPARSITY_PATTERN_H
#define __DOLFIN_SPARSITY_PATTERN_H

#include <dolfin/la/GenericSparsityPattern.h>

#include <dolfin/common/assert.h>
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
  SparsityPattern() = default;

  /// Create sparsity pattern for given global dimensions and local ranges.
  /// If range is a nullptr pointer the pattern is assumed to be serial.
  SparsityPattern(uint rank, uint const * dim, uint const * range = nullptr);

  /// Destructor
  ~SparsityPattern() override;

  //--- INTERFACE -------------------------------------------------------------

  /// Initialize with given tensor rank, global dimensions and local ranges.
  /// If range is a nullptr pointer the pattern is assumed to be serial
  void init(uint rank, uint const * dim, uint const * range = nullptr) override;

  /// Clear
  void clear() override;

  /// Insert non-zero entries
  void insert(uint const * num, uint const * const * idx) override;

  /// Return local size for given dimension
  auto size(uint i) const -> uint override;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  void apply() override;

  /// Is blocked
  auto is_blocked() const -> bool override;

  /// Is distributed
  auto is_distributed() const -> bool override;

  /// Return array with number of non-zeroes per local row
  void numNonZeroPerRow(uint nzrow[]) const override;

  /// Return array with number of non-zeroes per row for the given process rank
  /// and split between entries in the diagonal and off-diagonal portion of the
  /// matrix
  void numNonZeroPerRow(uint p_rank, uint d_nzrow[], uint o_nzrow[]) const override;

  /// Return total number of non-zeroes
  auto numNonZero() const -> uint override;

  /// Display sparsity pattern
  void disp() const override;

  //---------------------------------------------------------------------------

  /// Return array with row range for process rank
  void get_range(uint p_rank, uint range[]);

  /// Return number of local rows for process rank
  auto range_size(uint p_rank) const -> uint;

  ///
  void set_blocked();

private:

  /// Tensor rank
  uint rank_{0};

  /// Dimensions
  uint * dim_{nullptr};

  /// Range -array of size + 1 where size is size + 1:
  ///    range[rank], range[rank+1] is the range for processor
  uint ** range_{nullptr};

  /// Direct access to local range
  uint ** local_range_{nullptr};

  /// Flags
  bool initialized_{false};
  bool finalized_{false};
  bool blocked_{false};
  bool distributed_{false};

  /// Sparsity pattern represented as an array of sets.
  /// Each set corresponds to a row in the local range and contains the column
  /// positions of nonzero entries.

  /// Diagonal portion: submatrix such that row and column
  /// indices are in-range
  _ordered_set<uint> * d_entries_{nullptr};
  uint d_count_{0};

  /// Off-diagonal portion: entries such that only column indices are off-range
  _ordered_set<uint> * o_entries_{nullptr};
  uint o_count_{0};

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
inline auto SparsityPattern::is_blocked() const -> bool
{
  return blocked_;
}
//-----------------------------------------------------------------------------
inline auto SparsityPattern::is_distributed() const -> bool
{
  return distributed_;
}
//-----------------------------------------------------------------------------
inline auto SparsityPattern::size( uint i ) const -> uint
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
inline auto SparsityPattern::range_size( uint p_rank ) const -> uint
{
  dolfin_assert( distributed_ );
  return range_[0][p_rank + 1] - range_[0][p_rank];
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_SPARSITY_PATTERN_H */
