// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTIONAL_H
#define __DOLFIN_FUNCTIONAL_H

#include <dolfin/fem/Form.h>
#include <dolfin/fem/CoefficientMap.h>

namespace dolfin
{

class Functional : public Form
{

public:

  typedef Form::Coefficients Coefficients;

  static inline auto name() -> std::string { return "Functional"; }

  /// Constructor
  Functional(Mesh& mesh);

  /// Destructor
  ~Functional() override = default;

  /// Creator function
  template <class E> static inline
  auto create(Mesh& mesh, Coefficients& coefs) -> typename E::Functional *
  {
    return new typename E::Functional(mesh, coefs);
  }

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
