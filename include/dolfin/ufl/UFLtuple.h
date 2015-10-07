// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#ifndef __DOLFIN_UFL_TUPLE_H
#define __DOLFIN_UFL_TUPLE_H

#include <dolfin/ufl/UFLClass.h>

#include <iostream>
#include <sstream>
#include <string>

namespace ufl
{
class Expression;

/**
 *  DOCUMENTATION:
 *
 *  @class  tuple
 *
 *  @brief  Provides an interface complying with python tuple.
 */

template<class T>
class tuple : public Class
{

public:

  /// Constructor with default representation for given type
  tuple(T const& obj);

  /// Constructor with default representation for given type
  tuple(std::vector<T const *> const& objs);

  /// Constructor with default representation for given type
  tuple(tuple<T> const& other_tuple);

  ///
  tuple(repr_t const & repr);

  /// Create an empty tuple
  tuple();

  ///
  ~tuple();

  ///
  std::vector<Class const *> const operands(std::string const& name) const;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  ///
  dolfin::uint size() const;

  ///
  std::vector<T const *> const& operands() const;

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<T const *> const fill_objects(std::vector<repr_t> const& reprs);
  std::vector<T const *> const& objects() const;

  std::vector<T const *> const objects_;

  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_TUPLE_H */
