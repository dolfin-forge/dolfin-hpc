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
  GhostIterator & operator++();

  ///
  inline uint index() const;

  ///
  inline uint global_index() const;

  ///
  inline uint owner() const;

  ///
  inline bool valid() const;

  ///
  inline _set< uint > const & adj() const;

  ///
  template < class T >
  inline void adj_enqueue( Array< Array< T > > & container, T const & value ) const;

private:
  DistributedData const &            distdata_;
  _map< uint, uint >::const_iterator iter_;
};

//-----------------------------------------------------------------------------
inline GhostIterator & GhostIterator::operator++()
{
  ++iter_;
  return *this;
}
//-----------------------------------------------------------------------------
inline uint GhostIterator::index() const
{
  return iter_->first;
}
//-----------------------------------------------------------------------------
inline uint GhostIterator::global_index() const
{
  return distdata_.get_global( iter_->first );
}
//-----------------------------------------------------------------------------
inline uint GhostIterator::owner() const
{
  return iter_->second;
}
//-----------------------------------------------------------------------------
inline bool GhostIterator::valid() const
{
  return iter_ != distdata_.ghost_.end();
}
//-----------------------------------------------------------------------------
inline _set< uint > const & GhostIterator::adj() const
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
