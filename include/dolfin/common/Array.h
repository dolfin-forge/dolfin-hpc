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
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

#include <algorithm>
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
      std::vector<T>(),
      offset_(0),
      stride_(1)
  {
  }

  /// Create array of given size
  Array(uidx n) :
      std::vector<T>(n),
      offset_(0),
      stride_(1)
  {
  }

  /// Create array of given size with default value
  Array(uidx n, T const& t) :
      std::vector<T>(n, t),
      offset_(0),
      stride_(1)
  {
  }


  /// Create array given a range
  Array(T const * begin, T const * end) :
      std::vector<T>(begin, end),
      offset_(0),
      stride_(1)
  {
  }

  /// Create array given a range
  template<class Iterator>
  Array(Iterator const begin, Iterator const end) :
      std::vector<T>(begin, end),
      offset_(0),
      stride_(1)
  {
  }

  /// Copy constructor
  Array(Array<T> const& x) :
      std::vector<T>(x),
      offset_(x.offset_),
      stride_(x.stride_)
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

  /// Assignement operator
  Array<T>& operator=(Array<T> const& other)
  {
    if (this != &other)
    {
      std::vector<T>::operator=(other);
      offset_ = other.offset_;
      stride_ = other.stride_;
    }
    return *this;
  }

  /// Swap operator
  void swap(Array<T>& other)
  {
    if (this != &other)
    {
      std::vector<T>::swap(other);
      std::swap(offset_, other.offset_);
      std::swap(stride_, other.stride_);
    }
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

  /// Support missing assign() in SUN Studio
#ifdef _RWSTD_NO_MEMBER_TEMPLATES
  template<class Iterator>
  inline void assign(Iterator begin, Iterator end)
  {
    std::vector<T>::erase(this->begin(), this->end());
    std::copy(begin, end, std::back_inserter(*this));
  }
#endif


  /// Workaround defect in C++98: taking the address of the
  /// array of an empty std::vector is not allowed but it is often used in
  /// the code when dealing with MPI calls.
  /// In C++11 data() was added but the return value for an empty std::vector is
  /// left undefined by the standard. We do not want to rely on implementation
  /// details.
  inline T * ptr() { return (this->empty() ? NULL : &this->front()); }
  inline T const* ptr() const { return (this->empty() ? NULL : &this->front()); }

  /// Implement own semantics
  inline T * data() { return (this->empty() ? NULL : &this->front()); }
  inline T const* data() const { return (this->empty() ? NULL : &this->front()); }
  inline T * bound() { return this->data() + this->size(); }
  inline T const* bound() const { return this->data() + this->size(); }

  ///
  void operator%=(uint s)
  {
    if  (s == 0) { stride_ = this->size();  } else { stride_ = s; }
  }

  ///
  inline T * operator()(uint i)
  {
    dolfin_assert(i * stride_ < this->size());
    return &this->operator [](i * stride_);
  }

  ///
  inline uint& offset() { return offset_; }
  inline uint  offset() const { return offset_; }

  ///
  inline uint stride() const { return stride_; }

  ///
  inline uint dim(uint i) const
  {
      return (i == 0 ? this->size() / this->stride() :
              i == 1 ? this->stride() : 0);
  }

  /// Factor logic for array initialization
  inline static
  T * init(uidx n, T * src, T *& dst)
  {
    dolfin_assert(!(n == 0 && src != NULL));
    if (dst == NULL)
    {
      dst = (n > 0 ? new uint[n] : NULL);
    }
    if (src == NULL)
    {
      std::fill_n(dst, n, 0);
    }
    else
    {
      std::copy(src, src + n, dst);
    }
    return dst;
  }

  /// Factor logic for array initialization
  inline static
  T * init(uidx n, T * src)
  {
    dolfin_assert(!(n == 0 && src != NULL));
    T * dst = (n > 0 ? new uint[n] : NULL);
    if (src == NULL)
    {
      std::fill_n(dst, n, 0);
    }
    else
    {
      std::copy(src, src + n, dst);
    }
    return dst;
  }

  /// Dump data on the output
  void dump() const
  {
    for (typename Array<T>::const_iterator e = this->begin(); e != this->end();
         ++e)
    {
      cout << *e << "\n";
    }
  }

private:

  uint offset_;
  uint stride_;

};

//--- SPECIALIZATION ----------------------------------------------------------

template <class T>
class Array<T*> : public std::vector<T*>
{
public:

  /// Create empty array
  Array() :
      std::vector<T*>(),
      offset_(0),
      stride_(1)
  {
  }

  ///
  ~Array()
  {
  }

  /// Swap operator
  void swap(Array<T*>& other)
  {
    if (this != &other)
    {
      std::vector<T*>::swap(other);
      std::swap(offset_, other.offset_);
      std::swap(stride_, other.stride_);
    }
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

  ///
  void operator%=(uint s)
  {
    if  (s == 0) { stride_ = this->size();  } else { stride_ = s; }
  }

  ///
  inline uint& offset() { return offset_; }
  inline uint  offset() const { return offset_; }

  ///
  inline uint stride() const { return stride_; }

  ///
  inline uint dim(uint i) const
  {
      return (i == 0 ? this->size() / this->stride() :
              i == 1 ? this->stride() : 0);
  }

private:

  /// Disallow copy constructor
  Array(T const&) : offset_(0), stride_(0) {}

  /// Disallow assignement operator
  Array<T>& operator=(Array<T> const&) { return *this; }

  uint offset_;
  uint stride_;

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_ARRAY_H */
