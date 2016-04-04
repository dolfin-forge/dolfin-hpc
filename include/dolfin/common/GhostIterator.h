///
///
///

#ifndef __DOLFIN_GHOST_ITERATOR_H
#define __DOLFIN_GHOST_ITERATOR_H

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
  GhostIterator(DistributedData const& distdata) :
      distdata_(distdata),
      iter_(distdata_.ghosts_set_.begin())
  {
  }

  ///
  ~GhostIterator()
  {
  }

  ///
  GhostIterator& operator++()
  {
    ++iter_;
    return *this;
  }

  ///
  inline uint index() const
  {
    return *iter_;
  }

  ///
  inline uint global_index() const
  {
    return distdata_.get_global(*iter_);
  }

  ///
  inline uint owner() const
  {
    return distdata_.ghost_owner.find(*iter_)->second;
  }

  ///
  inline bool end() const
  {
    return iter_ == distdata_.ghosts_set_.end();
  }

  ///
  inline _set<uint> const& adj() const
  {
    return distdata_.shared_adj.find(*iter_)->second;
  }

private:

  DistributedData const& distdata_;
  _set<uint>::const_iterator iter_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_GHOST_ITERATOR_H */
