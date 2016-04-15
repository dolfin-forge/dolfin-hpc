// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2008-06-18
// Last changed: 2013-09-13

#include <dolfin/fem/BoundaryCondition.h>

#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type, Mesh& mesh,
                                     SubDomain const& sub_domain) :
    MeshDependent(mesh),
    type_(type),
    mesh_(mesh),
    sub_domain_index_(0),
    has_geometric_sub_domain_(true),
    geometric_sub_domain_(&sub_domain),
    sub_domain_markers_(NULL),
    local_sub_domain_markers_(true),
    sub_system_()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     MeshFunction<uint>& sub_domains,
                                     uint sub_domain) :
    MeshDependent(sub_domains.mesh()),
    type_(type),
    mesh_(sub_domains.mesh()),
    sub_domain_index_(sub_domain),
    has_geometric_sub_domain_(false),
    geometric_sub_domain_(NULL),
    sub_domain_markers_(&sub_domains),
    local_sub_domain_markers_(false),
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
    sub_domain_index_(0),
    has_geometric_sub_domain_(true),
    geometric_sub_domain_(&sub_domain),
    sub_domain_markers_(NULL),
    local_sub_domain_markers_(true),
    sub_system_(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     MeshFunction<uint>& sub_domains,
                                     uint sub_domain,
                                     SubSystem const sub_system) :
    MeshDependent(sub_domains.mesh()),
    type_(type),
    mesh_(sub_domains.mesh()),
    sub_domain_index_(sub_domain),
    has_geometric_sub_domain_(false),
    geometric_sub_domain_(NULL),
    sub_domain_markers_(&sub_domains),
    local_sub_domain_markers_(false),
    sub_system_(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::~BoundaryCondition()
{
  if (local_sub_domain_markers_)
  {
    delete sub_domain_markers_;
  }
}
//-----------------------------------------------------------------------------
void BoundaryCondition::init_markers(uint const entity_dim)
{
  // Create mesh function for sub domain markers on facets
  dolfin_assert(entity_dim < mesh().topology().dim());
  mesh().init(entity_dim);
  delete sub_domain_markers_;
  sub_domain_markers_ = new MeshFunction<uint>(mesh_, entity_dim);

  // Mark everything as sub domain 1
  (*sub_domain_markers_) = 1;

  // Mark the sub domain as sub domain 0
  geometric_sub_domain_->mark(*sub_domain_markers_, 0);
}
//-----------------------------------------------------------------------------
void BoundaryCondition::apply(GenericMatrix& A, GenericVector& b,
                              const BilinearForm& form,
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
                              const GenericVector& x, const BilinearForm& form,
                              SubSystem const sub_system)
{
  SubSystem defined = this->sub_system_;
  SubSystem enforced(sub_system, this->sub_system_);
  this->sub_system_ = enforced;
  apply(A, b, x, form);
  this->sub_system_ = defined;
}
//-----------------------------------------------------------------------------
std::string const& BoundaryCondition::type() const
{
  return type_;
}
//-----------------------------------------------------------------------------
Mesh& BoundaryCondition::mesh() const
{
  return mesh_;
}
//-----------------------------------------------------------------------------
bool BoundaryCondition::has_geometrical_sub_domain() const
{
  return has_geometric_sub_domain_;
}
//-----------------------------------------------------------------------------
SubDomain const& BoundaryCondition::sub_domain() const
{
  if (geometric_sub_domain_)
  {
    error("BoundaryCondition : no geometric subdomain defined");
  }
  return *geometric_sub_domain_;
}

//-----------------------------------------------------------------------------
uint const& BoundaryCondition::sub_domain_index() const
{
  return sub_domain_index_;
}
//-----------------------------------------------------------------------------
MeshFunction<uint> const& BoundaryCondition::sub_domain_markers() const
{
  if (sub_domain_markers_)
  {
    error("BoundaryCondition : no subdomain marker defined");
  }
  return *sub_domain_markers_;
}
//-----------------------------------------------------------------------------
SubSystem const& BoundaryCondition::sub_system() const
{
  return sub_system_;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
