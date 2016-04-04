///
///
///

#ifndef __DOLFIN_SHARED_ITERATOR_H
#define __DOLFIN_SHARED_ITERATOR_H

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
  SharedIterator(DistributedData const& distdata) :
      distdata_(distdata),
      iter_(distdata_.shared.begin())
  {
  }

  ///
  ~SharedIterator()
  {
  }

  ///
  SharedIterator& operator++()
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
  inline bool end() const
  {
    return iter_ == distdata_.shared.end();
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

#endif /* __DOLFIN_SHARED_ITERATOR_H */
