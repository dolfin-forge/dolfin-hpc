// Copyright (C) 2003 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2003-2007.
// Modified by Aurélien Larcher, 2014-217.
//
// First added:  2003-09-03
// Last changed: 2017-02-24

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

  /// Support for this construct does not exist in some STL implementations
  template<class Iterator>
  inline void append(Iterator begin, Iterator end)
  {
#ifdef __SUNPRO_CC
    for (Iterator it = begin; it != end; ++it) {this->push_back(*it); }
#else
    std::vector<T>::insert(this->end(), begin, end);
#endif
  }

};

//--- SPECIALIZATION ----------------------------------------------------------

template <class T>
class Array<T*> : public std::vector<T*>
{
public:

  /// Create empty array
  Array() :
      std::vector<T*>()
  {
  }

  ///
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

private:

  /// Disallow copy constructor
  Array(T const& other) {}

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_ARRAY_H */
