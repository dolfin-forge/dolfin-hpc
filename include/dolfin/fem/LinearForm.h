// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#ifndef __DOLFIN_LINEAR_FORM_H
#define __DOLFIN_LINEAR_FORM_H

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
  if (!test_space_)
  {
    test_space_ = this->create_space(0);
  }
  return *test_space_;
}

}

#endif
