// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETER_VALUE_H
#define __DOLFIN_PARAMETER_VALUE_H

#include <dolfin/common/types.h>

#include <string>

namespace dolfin
{

/// Base class for parameter values
struct ParameterValue
{

  /// Constructor
  ParameterValue();

  /// Destructor
  virtual ~ParameterValue();

  /// Equality
  bool operator==(ParameterValue const& other) const { return  this->cmp(other); }
  bool operator!=(ParameterValue const& other) const { return !this->cmp(other); }

  /// Assignment of int
  virtual ParameterValue const& operator=(int value);

  /// Assignment of uint
  virtual ParameterValue const& operator=(uint value);

  /// Assignment of real
  virtual ParameterValue const& operator=(real value);

  /// Assignment of bool
  virtual ParameterValue const& operator=(bool value);

  /// Assignment of string
  virtual ParameterValue const& operator=(std::string value);

  /// Cast to int
  virtual operator int() const;

  /// Cast to uint
  virtual operator uint() const;

  /// Cast to real
  virtual operator real() const;

  /// Cast to bool
  virtual operator bool() const;

  /// Cast to string
  virtual operator std::string() const;

  /// Value comparison
  virtual bool cmp(ParameterValue const& other) const = 0;

  /// Value comparison
  virtual ParameterValue const& operator>>(std::ostream& stream) const = 0;

  /// Name of value type
  virtual std::string type() const = 0;

  ///
  virtual ParameterValue * clone() const = 0;

};

///
template<class T>
struct parameter : public ParameterValue
{

  /// Constructor
  parameter(T value = T()) : ParameterValue(), v_(value) {}

  /// Copy constructor
  parameter(parameter<T> const& other) : ParameterValue(), v_(other.v_) {}

  /// Destructor
  ~parameter() {}

  /// Equality
  bool operator==(parameter const& other) const { return (v_ == other.v_); }
  bool operator!=(parameter const& other) const { return (v_ != other.v_); }

  /// Value comparison
  bool cmp(ParameterValue const& other) const { return v_ == static_cast<T>(other); }

  /// Output
  inline ParameterValue const& operator>>(std::ostream& stream) const
  {
    stream << v_; return *this;
  }

  /// Equality
  inline bool operator==(ParameterValue const& other) const
  {
    return (this->type() == other.type()) && (v_ == static_cast<T>(other));
  }
  bool operator!=(ParameterValue const& other) const { return !(*this == other); }

  /// Assignment
  ParameterValue& operator=(T value) { this->v_ = value; return *this; }

  /// Output
  template<class S>
  friend std::ostream& operator<<(std::ostream& stream, parameter<S> const& p);

  /// Cast
  operator T&()      { return v_; }
  operator T() const { return v_; }

  /// Name of value type
  std::string type() const;

  /// Clone
  ParameterValue * clone() const { return new parameter<T>(*this); }

private:

  T v_;

};

template<> inline std::string parameter<bool>::type()        const { return "bool"; }
template<> inline std::string parameter<int>::type()         const { return "int";  }
template<> inline std::string parameter<uint>::type()        const { return "uint"; }
template<> inline std::string parameter<real>::type()        const { return "real"; }
template<> inline std::string parameter<std::string>::type() const { return "string"; }

//-----------------------------------------------------------------------------
template<>
inline ParameterValue const& parameter<bool>::operator>>(std::ostream& stream) const
{
  stream << (v_ ? "true" : "false"); return *this;
}

//-----------------------------------------------------------------------------
template<>
inline ParameterValue const& parameter<std::string>::operator>>(std::ostream& stream) const
{
  stream << "\"" << v_ << "\""; return *this;
}

//-----------------------------------------------------------------------------
template<class T>
inline std::ostream& operator<<(std::ostream& stream, parameter<T> const& p)
{
  stream << p.v_; return stream;
}

//-----------------------------------------------------------------------------
template<class T>
inline std::ostream& operator<<(std::ostream& stream, ParameterValue const& p)
{
  p >> stream; return stream;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
