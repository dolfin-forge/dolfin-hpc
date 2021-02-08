// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BILINEAR_FORM_H
#define __DOLFIN_BILINEAR_FORM_H

#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class GenericMatrix;
class GenericVector;

class BilinearForm : public Form
{
public:
  using Coefficients = Form::Coefficients;

  static inline auto name() -> std::string
  {
    return "BilinearForm";
  }

  /// Constructor
  BilinearForm( Mesh & mesh );

  /// Destructor
  ~BilinearForm();

  /// Trial space
  auto trial_space() const -> FiniteElementSpace const &;

  /// Test space
  auto test_space() const -> FiniteElementSpace const &;

  /// Check whether linear system's dimensions match discrete spaces
  void check( GenericMatrix const & A, GenericVector const & b ) const;

  /// Creator function
  template < class E >
  static inline auto create( Mesh & mesh, Coefficients & coefs ) ->
    typename E::BilinearForm *
  {
    return new typename E::BilinearForm( mesh, coefs );
  }
};

//--- INLINES -----------------------------------------------------------------

inline auto BilinearForm::trial_space() const -> FiniteElementSpace const &
{
  return Form::spaces()[1];
}

//-----------------------------------------------------------------------------
inline auto BilinearForm::test_space() const -> FiniteElementSpace const &
{
  return Form::spaces()[0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
