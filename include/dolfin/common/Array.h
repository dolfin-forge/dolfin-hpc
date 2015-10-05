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

#include <iostream>

#include <dolfin/common/types.h>
#include <vector>

namespace dolfin
{

  /// Array is a container that provides O(1) access time to elements
  /// and O(1) memory overhead.
  ///
  /// It is a wrapper for std::vector, so see the STL manual for further
  /// details: http://www.sgi.com/tech/stl/

  template <class T>
  class Array : public std::vector<T>
  {
  public:

    /// Create empty array
    Array() : std::vector<T>() {}

    /// Create array of given size
    Array(uint n) : std::vector<T>(n) {}

    /// Create array of given size with default value
    Array(uint n, const T& t) : std::vector<T>(n, t) {}

    /// Create array containing two elements
    Array(const T& t0, const T& t1)
    {
      this->push_back(t0);
      this->push_back(t1);
    }

    /// Create array containing three elements
    Array(const T& t0, const T& t1, const T& t2)
    {
      this->push_back(t0);
      this->push_back(t1);
      this->push_back(t2);
    }

    /// Create array containing four elements
    Array(const T& t0, const T& t1, const T& t2, const T& t3)
    {
      this->push_back(t0);
      this->push_back(t1);
      this->push_back(t2);
      this->push_back(t3);
    }

    /// Create array containing five elements
    Array(const T& t0, const T& t1, const T& t2, const T& t3, const T& t4)
    {
      this->push_back(t0);
      this->push_back(t1);
      this->push_back(t2);
      this->push_back(t3);
      this->push_back(t4);
    }

    /// Copy constructor
    Array(const Array<T>& x) : std::vector<T>(x) {}

    /// Assign to all elements in the array
    const Array& operator=(const T& t)
    {
      for (uint i = 0; i < std::vector<T>::size(); ++i)
      (*this)[i] = t;
      return *this;
    }

    /// Destructor
    ~Array() {}

  };

}

#endif
