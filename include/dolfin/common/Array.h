// Copyright (C) 2003 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

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
  Array();

  /// Create array of given size
  Array(uidx n);

  /// Create array of given size with default value
  Array(uidx n, T const& t);


  /// Create array given a range
  Array(T const * begin, T const * end);

  /// Create array given a range
  template<class Iterator>
  Array(Iterator const begin, Iterator const end);

  /// Copy constructor
  Array(Array<T> const& x);

  /// Destructor
  ~Array();

  /// Assign to all elements in the array
  Array const& operator=(const T& t);

  /// Assignement operator
  Array<T>& operator=(Array<T> const& other);

  /// Swap operator
  void swap(Array<T>& other);

  /// Support for this construct does not exist in some STL implementations
  template<class Iterator>
  inline void append(Iterator begin, Iterator end);

  /// Support missing assign() in SUN Studio
#ifdef _RWSTD_NO_MEMBER_TEMPLATES
  template<class Iterator>
  inline void assign(Iterator begin, Iterator end);
#endif


  /// Workaround defect in C++98: taking the address of the
  /// array of an empty std::vector is not allowed but it is often used in
  /// the code when dealing with MPI calls.
  /// In C++11 data() was added but the return value for an empty std::vector is
  /// left undefined by the standard. We do not want to rely on implementation
  /// details.
  inline T * ptr();
  inline T const* ptr() const;

  /// Implement own semantics
  inline T * data();
  inline T const* data() const;

  inline T * bound();
  inline T const* bound() const;

  ///
  void operator%=(uint s);

  ///
  inline T * operator()(uint i);

  ///
  inline uint& offset();
  inline uint  offset() const;

  ///
  inline uint stride() const;

  ///
  inline uint dim(uint i) const;

  /// Factor logic for array initialization
  inline T * init(uidx n, T * src, T *& dst);

  /// Factor logic for array initialization
  inline T * init(uidx n, T * src);

  /// Dump data on the output
  void dump() const;

private:

  uint offset_;
  uint stride_;

};

/// Create empty array
template < class T >
Array< T >::Array()
  : std::vector< T >()
  , offset_( 0 )
  , stride_( 1 )
{
}

/// Create array of given size
template < class T >
Array< T >::Array( uidx n )
  : std::vector< T >( n )
  , offset_( 0 )
  , stride_( 1 )
{
}

/// Create array of given size with default value
template < class T >
Array< T >::Array( uidx n, T const & t )
  : std::vector< T >( n, t )
  , offset_( 0 )
  , stride_( 1 )
{
}

/// Create array given a range
template < class T >
Array< T >::Array( T const * begin, T const * end )
  : std::vector< T >( begin, end )
  , offset_( 0 )
  , stride_( 1 )
{
}

/// Create array given a range
template < class T >
template < class Iterator >
Array< T >::Array( Iterator const begin, Iterator const end )
  : std::vector< T >( begin, end )
  , offset_( 0 )
  , stride_( 1 )
{
}

/// Copy constructor
template < class T >
Array< T >::Array( Array< T > const & x )
  : std::vector< T >( x )
  , offset_( x.offset_ )
  , stride_( x.stride_ )
{
}

/// Destructor
template < class T >
Array< T >::~Array()
{
}

/// Assign to all elements in the array
template < class T >
Array< T > const & Array< T >::operator=( const T & t )
{
	std::fill( this->begin(), this->end(), t );
	return *this;
}

/// Assignement operator
template < class T >
Array< T > & Array< T >::operator=( Array< T > const & other )
{
	if ( this != &other )
	{
		std::vector< T >::operator=( other );

		offset_ = other.offset_;
		stride_ = other.stride_;
	}
	return *this;
}

template < class T >
void Array< T >::swap( Array< T > & other )
{
	if ( this != &other )
	{
		std::vector< T >::swap( other );
		std::swap( offset_, other.offset_ );
		std::swap( stride_, other.stride_ );
	}
}

/// Support for this construct does not exist in some STL implementations
template < class T >
template < class Iterator >
inline void Array< T >::append( Iterator begin, Iterator end )
{
#ifdef __SUNPRO_CC
	for ( Iterator it = begin; it != end; ++it )
	{
		this->push_back( *it );
	}
#else
	std::vector< T >::insert( this->end(), begin, end );
#endif
}

/// Support missing assign() in SUN Studio
#ifdef _RWSTD_NO_MEMBER_TEMPLATES
template < class T >
template < class Iterator >
inline void Array< T >::assign( Iterator begin, Iterator end )
{
	std::vector< T >::erase( this->begin(), this->end() );
	std::copy( begin, end, std::back_inserter( *this ) );
}
#endif

template < class T >
inline T * Array< T >::ptr()
{
	return ( this->empty() ? NULL : &this->front() );
}

template < class T >
inline T const * Array< T >::ptr() const
{
	return ( this->empty() ? NULL : &this->front() );
}

/// Implement own semantics
template < class T >
inline T * Array< T >::data()
{
	return ( this->empty() ? NULL : &this->front() );
}

template < class T >
inline T const * Array< T >::data() const
{
	return ( this->empty() ? NULL : &this->front() );
}

template < class T >
inline T * Array< T >::bound()
{
	return this->data() + this->size();
}

template < class T >
inline T const * Array< T >::bound() const
{
	return this->data() + this->size();
}

///
template < class T >
void Array< T >::operator%=( uint s )
{
	if ( s == 0 )
	{
		stride_ = this->size();
	}
	else
	{
		stride_ = s;
	}
}

///
template < class T >
inline T * Array< T >::operator()( uint i )
{
	dolfin_assert( i * stride_ < this->size() );
	return &this->operator[]( i * stride_ );
}

///
template < class T >
inline uint & Array< T >::offset()
{
	return offset_;
}

template < class T >
inline uint Array< T >::offset() const
{
	return offset_;
}

///
template < class T >
inline uint Array< T >::stride() const
{
	return stride_;
}

///
template < class T >
inline uint Array< T >::dim( uint i ) const
{
	return ( i == 0 ? this->size() / this->stride()
	                : i == 1 ? this->stride() : 0 );
}

/// Factor logic for array initialization
template < class T >
inline T * Array< T >::init( uidx n, T * src, T *& dst )
{
	dolfin_assert( !( n == 0 && src != NULL ) );
	if ( dst == NULL )
	{
		dst = ( n > 0 ? new uint[n] : NULL );
	}
	if ( src == NULL )
	{
		std::fill_n( dst, n, 0 );
	}
	else
	{
		std::copy( src, src + n, dst );
	}
	return dst;
}

/// Factor logic for array initialization
template < class T >
inline T * Array< T >::init( uidx n, T * src )
{
	dolfin_assert( !( n == 0 && src != NULL ) );
	T * dst = ( n > 0 ? new uint[n] : NULL );
	if ( src == NULL )
	{
		std::fill_n( dst, n, 0 );
	}
	else
	{
		std::copy( src, src + n, dst );
	}
	return dst;
}

/// Dump data on the output
template < class T >
void Array< T >::dump() const
{
	for ( typename Array< T >::const_iterator e = this->begin(); e != this->end();
	      ++e )
	{
		cout << *e << "\n";
	}
}

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
