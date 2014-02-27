// Copyright (C) 2007 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#include <dolfin/fem/BilinearForm.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BilinearForm::BilinearForm(Mesh& mesh) :
    Form(mesh),
    test_space_(NULL),
    trial_space_(NULL)
{}

//-----------------------------------------------------------------------------
BilinearForm::~BilinearForm()
{
  delete test_space_;
  delete trial_space_;
}

}
