// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_OWNED_ITERATOR_H
#define __DOLFIN_COMMON_OWNED_ITERATOR_H

#include <dolfin/common/DistributedData.h>

namespace dolfin
{

/**
 *  @class  OwnedIterator
 *
 *  @brief  Implements an iterator on owned entities for finalized distributed
 *          data only.
 */

class OwnedIterator
{

public:
  ///
  OwnedIterator( DistributedData const & distdata )
    : distdata_( distdata )
    , owner_( distdata.cached_ownership_ )
    , iter_( owner_.begin() )
  {
    if ( !distdata.is_finalized() )
    {
      error( "OwnedIterator : distributed data is not finalized" );
    }
  }

  ///
  ~OwnedIterator() = default;

  ///
  auto operator++() -> OwnedIterator &;

  ///
  inline auto index() const -> uint;

  ///
  inline auto global_index() const -> uint;

  ///
  inline auto is_shared() const -> uint;

  ///
  inline auto end() const -> bool;

private:
  DistributedData const &       distdata_;
  Array< uint > const &         owner_;
  Array< uint >::const_iterator iter_;
};

//-----------------------------------------------------------------------------
inline auto OwnedIterator::operator++() -> OwnedIterator &
{
  if ( iter_ == owner_.end() )
  {
    return *this;
  }
  ++iter_;
  while ( ( iter_ < owner_.end() ) && ( *iter_ != distdata_.pe_size_ )
          && ( *iter_ != distdata_.rank_ ) )
  {
    ++iter_;
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline auto OwnedIterator::index() const -> uint
{
  return iter_ - owner_.begin();
}
//-----------------------------------------------------------------------------
inline auto OwnedIterator::global_index() const -> uint
{
  return *iter_;
}
//-----------------------------------------------------------------------------
inline auto OwnedIterator::is_shared() const -> uint
{
  return ( *iter_ == distdata_.rank_ );
}
//-----------------------------------------------------------------------------
inline auto OwnedIterator::end() const -> bool
{
  return iter_ == owner_.end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_OWNED_ITERATOR_H */
