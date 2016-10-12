// Copyright (C) 2003 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2003-2007.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2003-09-03
// Last changed: 2007-04-24

#ifndef __DOLFIN_ARRAY_H
#define __DOLFIN_ARRAY_H

#include <dolfin/common/types.h>

#include <iostream>
#include <vector>

namespace dolfin
{

/// Array is a container that provides O(1) access time to elements
/// and O(1) memory overhead. => Thank you Captain Obvious!
///
/// It is a wrapper for std::vector, so see the STL manual for further
/// details: http://www.sgi.com/tech/stl/

template<class T>
class Array : public std::vector<T>
{
public:

  /// Create empty array
  Array() :
      std::vector<T>()
  {
  }

  /// Create array of given size
  Array(uint n) :
      std::vector<T>(n)
  {
  }

  /// Create array of given size with default value
  Array(uint n, T const& t) :
      std::vector<T>(n, t)
  {
  }

  /// Copy constructor
  Array(Array<T> const& x) :
      std::vector<T>(x)
  {
  }

  /// Assign to all elements in the array
  Array const& operator=(const T& t)
  {
    std::fill(this->begin(), this->end(), t);
    return *this;
  }

  /// Destructor
  ~Array()
  {
  }

  /// Cleanup array of allocated objects
  void free()
  {
    while (!this->empty())
    {
      delete this->back();
      this->pop_back();
    }
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_ARRAY_H */
