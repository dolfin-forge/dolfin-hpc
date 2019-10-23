// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_DICT_H
#define __DOLFIN_UFL_DICT_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLtype.h>

#include <iostream>
#include <sstream>
#include <string>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  tuple
 *
 *  @brief  Provides an interface complying with python dict.
 */

template<class KEY, class VALUE>
class dict : public Class
{

public:

  /// Constructor with default representation for given type
  dict(std::vector<std::pair<KEY const *, VALUE const *> > const& map);

  /// Constructor with default representation for given type
  dict(std::pair<KEY const *, VALUE const *> const& map);

  /// Constructor with default representation for given type
  dict(KEY const& key, VALUE const& value);

  /// Constructor with default representation for given type
  dict(dict const& other_dict);

  ///
  dict(repr_t const & repr);

  /// Create an empty dict
  dict();

  ///
  ~dict();

  ///
  std::vector<Class const *> const operands(std::string const& name) const;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  ///
  dolfin::uint size() const;

  ///
  std::vector<std::pair<KEY const *, VALUE const *> > const& map() const;

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  std::vector<std::pair<KEY const *, VALUE const *> > const fill_map(
      std::vector<repr_t> const& reprs);

  std::vector<std::pair<KEY const *, VALUE const *> > const map_;

  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_DICT_H */
