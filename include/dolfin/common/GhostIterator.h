// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_GHOST_ITERATOR_H
#define __DOLFIN_COMMON_GHOST_ITERATOR_H

#include <dolfin/common/DistributedData.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  @class  GhostIterator
 *
 *  @brief  Implements an iterator on ghost entities.
 */

class GhostIterator
{

public:
  ///
  GhostIterator( DistributedData const & distdata );

  ///
  ~GhostIterator() = default;

  ///
  auto operator++() -> GhostIterator &;

  ///
  inline auto index() const -> size_t;

  ///
  inline auto global_index() const -> size_t;

  ///
  inline auto owner() const -> size_t;

  ///
  inline auto valid() const -> bool;

  ///
  inline auto adj() const -> _set< size_t > const &;

  ///
  template < class T >
  inline auto adj_enqueue( std::vector< std::vector< T > > & container,
                           T const & value ) const -> void;

private:
  DistributedData const &                distdata_;
  _map< size_t, size_t >::const_iterator iter_;
};

//-----------------------------------------------------------------------------

inline GhostIterator::GhostIterator( DistributedData const & distdata )
  : distdata_( distdata )
  , iter_( distdata_.ghost_.begin() )
{
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::operator++() -> GhostIterator &
{
  ++iter_;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::index() const -> size_t
{
  return iter_->first;
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::global_index() const -> size_t
{
  return distdata_.get_global( iter_->first );
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::owner() const -> size_t
{
  return iter_->second;
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::valid() const -> bool
{
  return iter_ != distdata_.ghost_.end();
}

//-----------------------------------------------------------------------------

inline auto GhostIterator::adj() const -> _set< size_t > const &
{
  return distdata_.shared_.find( iter_->first )->second;
}

//-----------------------------------------------------------------------------

template < class T >
inline auto
  GhostIterator::adj_enqueue( std::vector< std::vector< T > > & container,
                              T const & value ) const -> void
{
  _set< size_t > const & a = distdata_.shared_.find( iter_->first )->second;
  for ( _set< size_t >::const_iterator it = a.begin(); it != a.end(); ++it )
  {
    container[*it].push_back( value );
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_GHOST_ITERATOR_H */
