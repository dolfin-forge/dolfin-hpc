// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETER_H
#define __DOLFIN_PARAMETER_H

#include <dolfin/common/types.h>

#include <string>

namespace dolfin
{

/// Base class for parameter values
class Parameter
{
public:
  /// Supported parameter types
  enum Type
  {
    bool_t   = 0,
    real_t   = 1,
    int_t    = 2,
    uint_t   = 3,
    string_t = 4
  };

public:
  /// Constructor
  Parameter( Type const & t )
    : type_( t )
  {
  }

  /// Destructor
  virtual ~Parameter()
  {
  }

  /// Return type of Parameter
  Type type() const
  {
    return type_;
  }

  ///
  // virtual Parameter * clone() const = 0;
protected:
  Type type_;
};

/// specialization class which actually holds the data
template < class T >
struct parameter : public Parameter
{
  /// Constructor
  parameter( T const & value, Type const & t )
    : Parameter( t )
    , v_( value )
  {
  }

  /// Copy constructor
  parameter( parameter< T > const & other )
    : Parameter( other.type_ )
    , v_( other.v_ )
  {
  }

  /// Destructor
  ~parameter()
  {
  }

  /// Assignment
  Parameter & set( T const & value )
  {
    this->v_ = value;
    return *this;
  }

  /// Access
  T & get()
  {
    return v_;
  }

  T const & get() const
  {
    return v_;
  }

  /// Clone
  // Parameter * clone() const { return new parameter<T>(*this); }

private:
  T v_;
};
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
