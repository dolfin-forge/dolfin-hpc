// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_LIST_H
#define __DOLFIN_UFL_LIST_H

#include <dolfin/ufl/UFLIntegral.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  List
 *
 *  @brief  Provides an interface complying with Lists.
 */

class List : public Class
{
public:

  ///
  List(std::vector<Integral const *> const& integrals);

  ///
  List(repr_t const & repr);

  ///
  ~List() override;

  std::vector<Integral const *> const& get_integrals() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  std::vector<Integral const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Integral const *> integrals_;

  repr_t const repr_;
  std::string const str_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_LIST_H */
