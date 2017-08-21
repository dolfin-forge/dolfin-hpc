// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-02-22
//

#ifndef __DOLFIN_CLONABLE_H
#define __DOLFIN_CLONABLE_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace dolfin
{

template <class T>
struct Clonable
{

  /// Clone instance
  inline T * clone() const { return new T(static_cast<T const&>(*this)); };

  /// Clone second argument into first argument if non-NULL
  static void clone(Clonable<T> *& p0, Clonable<T> const * p1)
  {
    dolfin_assert(p0 == NULL);
    if (p1 != NULL) { p0 = p1->clone(); }
  }

};

/// Clone second argument into first argument if non-NULL
template <class CloneT>
static void cloneptr(CloneT *& p0, CloneT const * p1)
{
  dolfin_assert(p0 == NULL);
  if (p1 != NULL) { p0 = p1->clone(); }
}

/// Clone second argument into first argument if non-NULL
template <class CloneT>
static void cloneptr(CloneT const *& p0, CloneT const * p1)
{
  dolfin_assert(p0 == NULL);
  if (p1 != NULL) { p0 = p1->clone(); }
}

} /* namespace dolfin */

#endif /* __DOLFIN_CLONABLE_H */
