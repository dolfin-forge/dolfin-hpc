// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_OBJECT_H_
#define __UFL_OBJECT_H_

#include <string>
#include <vector>

#include <dolfin/ufl/UFLrepr.h>

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

  /// __eq__
  virtual bool operator == (Object const& other) const;

protected:

  ///
  Object() {}

  ///
  virtual ~Object() {}

};

//-----------------------------------------------------------------------------
inline bool Object::operator == (Object const& other) const
{
    return ( other.repr() == this->repr() );
}

} /* namespace ufl */
#endif /* __UFL_OBJECT_H_ */
