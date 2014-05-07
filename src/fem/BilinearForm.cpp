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
{
}

//-----------------------------------------------------------------------------
BilinearForm::~BilinearForm()
{
  delete test_space_;
  delete trial_space_;
}

//-----------------------------------------------------------------------------
void BilinearForm::check(GenericMatrix const& A, GenericVector const& b) const
{
  uint const M = this->test_space().dofmap().global_dimension();
  uint const N = this->trial_space().dofmap().global_dimension();

  if (M != A.size(0))
    error("Incorrect dimension 0 of matrix for given test space.");
  if (N != A.size(1))
    error("Incorrect dimension 1 of matrix for given trial space");
  if (M != b.size())
    error("Incorrect dimension of vector for given test space");
}

}
