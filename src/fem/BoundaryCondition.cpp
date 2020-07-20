// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/BoundaryCondition.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type, Mesh& mesh,
                                     SubDomain const& sub_domain) :
    MeshDependent(mesh),
    type_(type),
    mesh_(mesh),
    sub_domain_(sub_domain),
    sub_system_()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type, Mesh& mesh,
                                     SubDomain const& sub_domain,
                                     SubSystem const sub_system) :
    MeshDependent(mesh),
    type_(type),
    mesh_(mesh),
    sub_domain_(sub_domain),
    sub_system_(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::~BoundaryCondition() = default;
//-----------------------------------------------------------------------------
void BoundaryCondition::apply(GenericMatrix& A, GenericVector& b,
                              BilinearForm const& form,
                              SubSystem const sub_system)
{
  SubSystem defined = this->sub_system_;
  SubSystem enforced(sub_system, this->sub_system_);
  this->sub_system_ = enforced;
  apply(A, b, form);
  this->sub_system_ = defined;
}
//-----------------------------------------------------------------------------
void BoundaryCondition::apply(GenericMatrix& A, GenericVector& b,
                              GenericVector const& x, BilinearForm const& form,
                              SubSystem const sub_system)
{
  SubSystem defined = this->sub_system_;
  SubSystem enforced(sub_system, this->sub_system_);
  this->sub_system_ = enforced;
  apply(A, b, x, form);
  this->sub_system_ = defined;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
