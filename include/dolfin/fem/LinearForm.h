// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LINEAR_FORM_H
#define __DOLFIN_LINEAR_FORM_H

#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>

namespace dolfin
{

class LinearForm : public Form
{

public:
  typedef Form::Coefficients Coefficients;

  static inline auto name() -> std::string
  {
    return "LinearForm";
  }

  /// Constructor
  LinearForm( Mesh & mesh );

  /// Destructor
  ~LinearForm();

  /// Test space
  auto test_space() const -> FiniteElementSpace const &;

  /// Creator function
  template < class E >
  static inline auto create( Mesh & mesh, CoefficientMap & coefs ) ->
    typename E::LinearForm *
  {
    return new typename E::LinearForm( mesh, coefs );
  }
};

//--- INLINES -----------------------------------------------------------------

inline auto LinearForm::test_space() const -> FiniteElementSpace const &
{
  return Form::spaces()[0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
