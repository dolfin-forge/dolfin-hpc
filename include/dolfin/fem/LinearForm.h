// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LINEAR_FORM_H
#define __DOLFIN_LINEAR_FORM_H

#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/CoefficientMap.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class LinearForm : public Form
{

public:

  typedef Form::Coefficients Coefficients;

  static inline std::string name() { return "LinearForm"; }

  /// Constructor
  LinearForm(Mesh& mesh);

  /// Destructor
  ~LinearForm() override;

  /// Test space
  FiniteElementSpace const& test_space() const;

  /// Creator function
  template <class E> static inline
  typename E::LinearForm * create(Mesh& mesh, CoefficientMap& coefs)
  {
    return new typename E::LinearForm(mesh, coefs);
  }

private:

  mutable FiniteElementSpace * test_space_;

};

//--- INLINES -----------------------------------------------------------------

inline FiniteElementSpace const& LinearForm::test_space() const
{
  if (!test_space_)
  {
    test_space_ = this->create_space(0);
  }
  return *test_space_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif
