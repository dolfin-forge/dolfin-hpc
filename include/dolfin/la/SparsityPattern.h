// Copyright (C) 2007 Garth N. Wells
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SPARSITY_PATTERN_H
#define __DOLFIN_SPARSITY_PATTERN_H

#include <dolfin/la/GenericSparsityPattern.h>

#include <dolfin/common/assert.h>
#include <dolfin/common/types.h>

#include <vector>

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
  SparsityPattern(size_t rank, size_t const * dim, size_t const * range = nullptr);

  /// Destructor
  ~SparsityPattern() override;

  //--- INTERFACE -------------------------------------------------------------

  /// Initialize with given tensor rank, global dimensions and local ranges.
  /// If range is a nullptr pointer the pattern is assumed to be serial
  void init(size_t rank, size_t const * dim, size_t const * range = nullptr) override;

  /// Clear
  void clear() override;

  /// Insert non-zero entries
  void insert(size_t const * num, size_t const * const * idx) override;

  /// Return local size for given dimension
  auto size(size_t i) const -> size_t override;

  /// Finalize sparsity pattern (needed by most parallel la backends)
  void apply() override;

  /// Is blocked
  auto is_blocked() const -> bool override;

  /// Is distributed
  auto is_distributed() const -> bool override;

  /// Return array with number of non-zeroes per local row
  void numNonZeroPerRow(size_t nzrow[]) const override;

  /// Return array with number of non-zeroes per row for the given process rank
  /// and split between entries in the diagonal and off-diagonal portion of the
  /// matrix
  void numNonZeroPerRow(size_t p_rank, size_t d_nzrow[], size_t o_nzrow[]) const override;

  /// Return total number of non-zeroes
  auto numNonZero() const -> size_t override;

  /// Display sparsity pattern
  void disp() const override;

  //---------------------------------------------------------------------------

  /// Return array with row range for process rank
  auto get_range(size_t p_rank, size_t range[]) const -> void;

  /// Return number of local rows for process rank
  auto range_size(size_t p_rank) const -> size_t;

  ///
  auto set_blocked() -> void;

  //---------------------------------------------------------------------------

  // access diagonal columns of the matrix for a given row
  auto diagonal_entries( size_t row, std::vector< size_t > & cols ) const -> void;

  // access off-diagonal columns of the matrix for a given row
  auto off_diagonal_entries( size_t row, std::vector< size_t > & cols ) const -> void;

private:
  size_t const pe_size;
  size_t const pe_rank;

  /// Tensor rank
  size_t rank_{0};

  /// Dimensions
  std::vector< size_t > dim_;

  /// Range -array of size + 1 where size is size + 1:
  ///    range[rank], range[rank+1] is the range for processor
  std::vector< std::vector< size_t > > range_;

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
  std::vector< _ordered_set<size_t> > d_entries_;
  size_t d_count_{0};

  /// Off-diagonal portion: entries such that only column indices are off-range
  std::vector< _ordered_set<size_t> > o_entries_;
  size_t o_count_{0};

  /// Additionally provide data structure to store remote entries i,e such that
  /// row indices are not in-range
  _ordered_map<size_t, _ordered_set<size_t> > r_entries_;

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

inline auto SparsityPattern::size( size_t i ) const -> size_t
{
  return ( range_[i][pe_rank+1] - range_[i][pe_rank] );
}

//-----------------------------------------------------------------------------

inline void SparsityPattern::get_range( size_t p_rank, size_t range[] ) const
{
  dolfin_assert( p_rank < pe_size );
  // For a serial pattern p_rank is only zero
  range[0] = range_[0][p_rank];
  range[1] = range_[0][p_rank + 1];
}

//-----------------------------------------------------------------------------

inline auto SparsityPattern::range_size( size_t p_rank ) const -> size_t
{
  dolfin_assert( p_rank < pe_size );
  return range_[0][p_rank + 1] - range_[0][p_rank];
}

//-----------------------------------------------------------------------------

inline auto SparsityPattern::diagonal_entries( size_t row, std::vector< size_t > & cols ) const -> void
{
  dolfin_assert( row < this->size(0) );
  cols.resize( d_entries_[row].size() );
  std::copy( d_entries_[row].begin(), d_entries_[row].end(), cols.begin() );
}

//-----------------------------------------------------------------------------

inline auto SparsityPattern::off_diagonal_entries( size_t row, std::vector< size_t > & cols ) const -> void
{
  dolfin_assert( row < this->size(0) );
  cols.resize( d_entries_[row].size() );
  std::copy( d_entries_[row].begin(), d_entries_[row].end(), cols.begin() );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_SPARSITY_PATTERN_H */
