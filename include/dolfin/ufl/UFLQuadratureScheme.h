// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_QUADRATURE_SCHEME_H_
#define __UFL_QUADRATURE_SCHEME_H_

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
  ~QuadratureScheme();

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  std::string const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __UFL_QUADRATURE_SCHEME_H_ */
