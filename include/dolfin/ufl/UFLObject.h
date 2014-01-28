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
    std::vector<Object const *> const& prototype) const
{
  std::stringstream ret;
  std::vector<Object const *>::const_iterator arg = prototype.begin();
  ret << (*arg)->repr();
  for ( ++arg; arg != prototype.end(); ++arg)
  {
    ret << ", " << (*arg)->repr();
  }
  return ret.str();
}

} /* namespace ufl */
#endif /* __UFL_OBJECT_H_ */
