// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CLASS_H_
#define __UFL_CLASS_H_

#include <string>
#include <vector>

#include <dolfin/ufl/UFLObject.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Class
 *
 *  @brief  Provides an interface for Python classes from UFL.
 */

class Class : public Object
{

public:

  ///
  std::string const& name() const;

  /// __repr__
  virtual repr_t const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

  /// __eq__
  virtual bool operator == (Class const& other) const;

protected:

  Class();
  Class(std::string const& name);

  virtual ~Class();

  ///
  virtual void display() const;

private:

  std::string const name_;
  static repr_t const default_repr_;
  static std::string const default_str_;

};

//-----------------------------------------------------------------------------
inline bool Class::operator == (Class const& other) const
{
    return ( other.repr() == this->repr() );
}

//-----------------------------------------------------------------------------
class ValueArray : public std::vector<uint>
{

public:

  ///
  ValueArray();

  ///
  ValueArray(uint const i);

  ///
  ValueArray(uint const k, uint const i);

  ///
  ~ValueArray();

  ///
  std::string const str() const;

private:

};


} /* namespace ufl */
#endif /* __UFL_CLASS_H */
