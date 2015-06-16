// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#ifndef __UFL_ARRAY_H_
#define __UFL_ARRAY_H_

#include <dolfin/ufl/UFLObject.h>

#include <iostream>
#include <sstream>
#include <string>

namespace ufl
{

class array : public Object
{

public:

  /// Constructor with default representation for given type
  array(Object const& other, bool const unpack = false);

  ///
  ~array();

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

  ///
  repr_t make_repr(std::vector<Object const *> const& prototype) const;

  ///
  static Object::repr_t unpack(Object::repr_t const& repr_array);

protected:

private:

  Object const& obj_;
  bool const unpack_;
  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace icorne */
#endif /* __UFL_ARRAY_H_ */
