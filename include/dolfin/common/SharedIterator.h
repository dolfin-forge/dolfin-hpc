// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_SHARED_ITERATOR_H
#define __DOLFIN_COMMON_SHARED_ITERATOR_H

#include <dolfin/common/DistributedData.h>

namespace dolfin
{

/**
 *  @class  SharedIterator
 *
 *  @brief  Implements an iterator on shared entities.
 */

class SharedIterator
{

public:
  ///
  SharedIterator( DistributedData const & distdata )
    : distdata_( distdata )
    , iter_( distdata_.shared_.begin() )
  {
  }

  ///
  ~SharedIterator() = default;

  ///
  auto operator++() -> SharedIterator &;

  ///
  inline auto index() const -> uint;

  ///
  inline auto global_index() const -> uint;

  ///
  inline auto owner() const -> uint;

  ///
  inline auto is_owned() const -> bool;

  ///
  inline auto valid() const -> bool;

  ///
  inline auto adj() const -> _set< uint > const &;

  ///
  template < typename T >
  inline void adj_enqueue( Array< Array< T > > & container, T const & value ) const;

private:
  DistributedData const &                    distdata_;
  DistributedData::SharedSet::const_iterator iter_;
};

//-----------------------------------------------------------------------------
inline auto SharedIterator::operator++() -> SharedIterator &
{
  ++iter_;
  return *this;
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::index() const -> uint
{
  return iter_->first;
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::global_index() const -> uint
{
  return distdata_.get_global( iter_->first );
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::owner() const -> uint
{
  return distdata_.get_owner( iter_->first );
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::is_owned() const -> bool
{
  return distdata_.is_owned( iter_->first );
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::valid() const -> bool
{
  return iter_ != distdata_.shared_.end();
}
//-----------------------------------------------------------------------------
inline auto SharedIterator::adj() const -> _set< uint > const &
{
  return iter_->second;
}
//-----------------------------------------------------------------------------
template < typename T >
inline void SharedIterator::adj_enqueue( Array< Array< T > > & container,
                                         T const & value ) const
{
  _set< uint > const & a = iter_->second;
  for ( _set< uint >::const_iterator it = a.begin(); it != a.end(); ++it )
  {
    container[*it].push_back( value );
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_SHARED_ITERATOR_H */
