// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_OBJECT_H_
#define __UFL_OBJECT_H_

#include <dolfin/ufl/UFLrepr.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Object
 *
 *  @brief  Provides an interface for Python objects from UFL.
 */

class Object
{

public:

  typedef ufl::repr repr_t;

  /// __repr__
  virtual repr_t const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

  ///
  virtual repr_t const make_repr(
      std::vector<Object const *> const& prototype) const = 0;

  /// __eq__
  virtual bool operator ==(Object const& other) const;

protected:

  ///
  Object()
  {
  }

  ///
  virtual ~Object()
  {
  }

  ///
  virtual void display() const = 0;

  ///
  virtual std::vector<repr_t> const make_args_repr(
      repr_t const& repr) const = 0;

  ///
  std::vector<Object const *> make_args(std::vector<repr_t> const& repr) const;

  /// Create from representation
  static Object * create(repr_t representation);

};

//-----------------------------------------------------------------------------
inline bool Object::operator ==(Object const& other) const
{
  return (other.repr() == this->repr());
}

//-----------------------------------------------------------------------------
inline void Object::display() const
{
  std::cout << "Object" << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//-----------------------------------------------------------------------------
inline Object::repr_t const Object::make_repr(
    std::vector<Object const *> const& args) const
{
  std::stringstream ret;
  std::vector<Object const *>::const_iterator arg = args.begin();
  ret << (*arg)->repr();
  for (++arg; arg != args.end(); ++arg)
  {
    ret << ", " << (*arg)->repr();
  }
  return ret.str();
}

//-----------------------------------------------------------------------------
inline std::vector<Object const *> Object::make_args(
    std::vector<repr_t> const& repr) const
{
  std::vector<Object const *> args;
  for (std::vector<repr_t>::const_iterator it = repr.begin(); it != repr.end();
      ++it)
  {
    args.push_back(this);
  }
  return args;
}

//-----------------------------------------------------------------------------
inline std::vector<Object::repr_t> const Object::make_args_repr(
    repr_t const& repr) const
{
  std::vector<Object::repr_t> args;
  std::string str = repr;
  std::string delimiter = ", ";

  size_t scpos = 0;
  size_t currpos = 0;
  size_t openbrace = 0;
  std::string token;
  while (scpos != std::string::npos)
  {
    scpos = str.find(delimiter, currpos);
    openbrace = str.find("(", currpos);
    // Take into account of
    // - one argument
    // - class arguments
    // - tuple argument
    // - dict argument
    if (scpos == std::string::npos || openbrace > scpos
        || (openbrace < scpos && str.find(")", currpos) < scpos)
        || (str.find("{", currpos) < scpos && str.find("}", currpos) < scpos))
    {
      token = str.substr(0, scpos);
      args.push_back(Object::repr_t(token));
      str.erase(0, scpos + delimiter.length());
      currpos = 0;
    }
    else
    {
      currpos = scpos + 1;
    }
  }
  return args;
}

} /* namespace ufl */
#endif /* __UFL_OBJECT_H_ */
