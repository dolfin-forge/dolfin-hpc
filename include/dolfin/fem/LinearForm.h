// Copyright (C) 2007 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#ifndef __LINEAR_FORM_H
#define __LINEAR_FORM_H

#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElementSpace.h>

#include <ufc.h>

namespace dolfin
{

class LinearForm : public Form
{
public:

  /// Constructor
  LinearForm(Mesh& mesh);

  /// Destructor
  ~LinearForm();

  /// Test space
  FiniteElementSpace const& test_space() const;

private:

  mutable FiniteElementSpace * test_space_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline FiniteElementSpace const& LinearForm::test_space() const
{
  if(!test_space_)
  {
    ufc::finite_element * test = this->form().create_finite_element(0);
    test_space_ = new FiniteElementSpace(mesh(), *test, true);
  }
  return *test_space_;
}

}

#endif
