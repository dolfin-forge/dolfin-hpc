// Copyright (C) 2007 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/BilinearForm.h>

#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BilinearForm::BilinearForm(Mesh& mesh) :
    Form(mesh),
    test_space_(nullptr),
    trial_space_(nullptr)
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
  // FIXME num_entities should maybe be stored somewhere else
  std::vector< size_t > num_entities_M(
    this->test_space().element()->topological_dimension(), 0 );
  std::vector< size_t > num_entities_N(
    this->trial_space().element()->topological_dimension(), 0 );

  for ( uint d = 0; d <= num_entities_M.size(); ++d )
    if ( this->mesh().topology().connectivity( d ) )
      num_entities_M[d] = this->mesh().topology().global_size( d );

  for ( uint d = 0; d <= num_entities_N.size(); ++d )
    if ( this->mesh().topology().connectivity( d ) )
      num_entities_N[d] = this->mesh().topology().global_size( d );

  uint const M = this->test_space().dofmap().global_dimension( num_entities_M );
  uint const N = this->trial_space().dofmap().global_dimension( num_entities_N );

  if(A.size(0) != A.size(1))
  {
    error("Only square linear system are supported.");
  }
  if (M != A.size(0))
  {
    error("Incorrect dimension 0 of matrix for given test space.\n"
          "Space  : %u ; Matrix : %u", M, A.size(0));
  }
  if (N != A.size(1))
  {
    error("Incorrect dimension 1 of matrix for given trial space\n"
          "Space  : %u ; Matrix : %u", N, A.size(1));
  }
  if (M != b.size())
  {
    error("Incorrect dimension of vector for given test space\n"
          "Space  : %u ; Vector : %u", M, b.size(0));
  }
}

} /* namespace dolfin */
