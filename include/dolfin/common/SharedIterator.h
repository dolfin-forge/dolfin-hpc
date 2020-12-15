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
  SharedIterator & operator++();

  ///
  inline uint index() const;

  ///
  inline uint global_index() const;

  ///
  inline uint owner() const;

  ///
  inline bool is_owned() const;

  ///
  inline bool valid() const;

  ///
  inline _set< uint > const & adj() const;

  ///
  template < typename T >
  inline void adj_enqueue( Array< Array< T > > & container, T const & value ) const;

private:
  DistributedData const &                    distdata_;
  DistributedData::SharedSet::const_iterator iter_;
};

//-----------------------------------------------------------------------------
inline SharedIterator & SharedIterator::operator++()
{
  ++iter_;
  return *this;
}
//-----------------------------------------------------------------------------
inline uint SharedIterator::index() const
{
  return iter_->first;
}
//-----------------------------------------------------------------------------
inline uint SharedIterator::global_index() const
{
  return distdata_.get_global( iter_->first );
}
//-----------------------------------------------------------------------------
inline uint SharedIterator::owner() const
{
  return distdata_.get_owner( iter_->first );
}
//-----------------------------------------------------------------------------
inline bool SharedIterator::is_owned() const
{
  return distdata_.is_owned( iter_->first );
}
//-----------------------------------------------------------------------------
inline bool SharedIterator::valid() const
{
  return iter_ != distdata_.shared_.end();
}
//-----------------------------------------------------------------------------
inline _set< uint > const & SharedIterator::adj() const
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
