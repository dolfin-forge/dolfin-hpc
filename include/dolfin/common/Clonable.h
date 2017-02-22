// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-02-22
//

#ifndef __DOLFIN_CLONABLE_H
#define __DOLFIN_CLONABLE_H

namespace dolfin
{

template <class T>
struct Clonable
{

  inline T * clone() const { return new T(static_cast<T const&>(*this)); };

};

} /* namespace dolfin */

#endif /* __DOLFIN_CLONABLE_H */
