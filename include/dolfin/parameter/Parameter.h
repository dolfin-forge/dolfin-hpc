// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2005-03-24

#ifndef __DOLFIN_PARAMETER_H
#define __DOLFIN_PARAMETER_H

#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

class ParameterValue;

/// This class represents a parameter of some given type.
/// Supported value types are bool. int, uint, real, and string.

class Parameter
{
public:

  /// Supported parameter types
  enum Type
  {
    type_real, type_int, type_uint, type_bool, type_string
  };

  /// Create int-valued parameter
  Parameter(int value);

  /// Create int-valued parameter
  Parameter(uint value);

  /// Create real-valued parameter
  Parameter(real value);

  /// Create bool-valued parameter
  Parameter(bool value);

  /// Create string-valued parameter
  Parameter(std::string value);

  /// Create string-valued parameter
  Parameter(const char* value);

  /// Copy constructor
  Parameter(Parameter const& parameter);

  /// Destructor
  ~Parameter();

  /// Equality
  bool operator==(Parameter const& other) const;
  bool operator!=(Parameter const& other) const;

  /// Assignment of int
  Parameter const& operator=(int value);

  /// Assignment of uint
  Parameter const& operator=(uint value);

  /// Assignment of real
  Parameter const& operator=(real value);

  /// Assignment of bool
  Parameter const& operator=(bool value);

  /// Assignment of string
  Parameter const& operator=(std::string value);

  /// Assignment of Parameter
  Parameter const& operator=(Parameter const& parameter);

  /// Cast parameter to int
  operator int() const;

  /// Cast parameter to uint
  operator uint() const;

  /// Cast parameter to real
  operator real() const;

  /// Cast parameter to bool
  operator bool() const;

  /// Cast parameter to string
  operator std::string() const;

  /// Return type of parameter
  Type type() const;

  /// Output
  friend LogStream& operator<<(LogStream& stream, Parameter const& parameter);

  /// Output
  friend std::ostream& operator<<(std::ostream& stream, Parameter const& parameter);

  /// Friends
  friend class XMLFile;

private:

  // Pointer to parameter value
  ParameterValue* value_;

  // Type of parameter
  Type type_;

};

LogStream& operator<<(LogStream& stream, Parameter const& parameter);

}

#endif
