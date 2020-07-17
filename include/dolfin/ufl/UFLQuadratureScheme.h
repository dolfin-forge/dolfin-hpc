// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_QUADRATURE_SCHEME_H
#define __DOLFIN_UFL_QUADRATURE_SCHEME_H

#include <dolfin/ufl/UFLClass.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLQuadratureScheme
 *
 *  @brief  Provides an interface for UFL QuadratureScheme.
 */

class QuadratureScheme : public Class
{

public:

  ///
  QuadratureScheme();

  ///
  ~QuadratureScheme() override;

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

private:

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_QUADRATURE_SCHEME_H */
