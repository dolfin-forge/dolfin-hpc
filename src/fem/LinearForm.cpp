// Copyright (C) 2007 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#include <dolfin/fem/LinearForm.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
LinearForm::LinearForm(Mesh& mesh) :
    Form(mesh),
    test_space_(NULL)
{}

//-----------------------------------------------------------------------------
LinearForm::~LinearForm()
{
  delete test_space_;
}

}
