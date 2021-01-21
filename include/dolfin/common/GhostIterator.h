// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_GHOST_ITERATOR_H
#define __DOLFIN_COMMON_GHOST_ITERATOR_H

#include <dolfin/common/DistributedData.h>

namespace dolfin
{

/**
 *  @class  GhostIterator
 *
 *  @brief  Implements an iterator on ghost entities.
 */

class GhostIterator
{

public:
  ///
  GhostIterator( DistributedData const & distdata )
    : distdata_( distdata )
    , iter_( distdata_.ghost_.begin() )
  {
  }

  ///
  ~GhostIterator() = default;

  ///
  auto operator++() -> GhostIterator &;

  ///
  inline auto index() const -> uint;

  ///
  inline auto global_index() const -> uint;

  ///
  inline auto owner() const -> uint;

  ///
  inline auto valid() const -> bool;

  ///
  inline auto adj() const -> _set< uint > const &;

  ///
  template < class T >
  inline void adj_enqueue( Array< Array< T > > & container, T const & value ) const;

private:
  DistributedData const &            distdata_;
  _map< uint, uint >::const_iterator iter_;
};

//-----------------------------------------------------------------------------
inline auto GhostIterator::operator++() -> GhostIterator &
{
  ++iter_;
  return *this;
}
//-----------------------------------------------------------------------------
inline auto GhostIterator::index() const -> uint
{
  return iter_->first;
}
//-----------------------------------------------------------------------------
inline auto GhostIterator::global_index() const -> uint
{
  return distdata_.get_global( iter_->first );
}
//-----------------------------------------------------------------------------
inline auto GhostIterator::owner() const -> uint
{
  return iter_->second;
}
//-----------------------------------------------------------------------------
inline auto GhostIterator::valid() const -> bool
{
  return iter_ != distdata_.ghost_.end();
}
//-----------------------------------------------------------------------------
inline auto GhostIterator::adj() const -> _set< uint > const &
{
  return distdata_.shared_.find( iter_->first )->second;
}
//-----------------------------------------------------------------------------
template < class T >
inline void GhostIterator::adj_enqueue( Array< Array< T > > & container,
                                        T const & value ) const
{
  _set< uint > const & a = distdata_.shared_.find( iter_->first )->second;
  for ( _set< uint >::const_iterator it = a.begin(); it != a.end(); ++it )
  {
    container[*it].push_back( value );
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_GHOST_ITERATOR_H */
