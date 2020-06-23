// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This abstraction was rewritten to violate abstraction the least possible
// instead of "in every possible way".
// Additionally support for subsystems was added.

#ifndef __DOLFIN_BOUNDARY_CONDITION_H
#define __DOLFIN_BOUNDARY_CONDITION_H

#include <dolfin/mesh/MeshDependent.h>

#include <dolfin/evolution/Time.h>
#include <dolfin/fem/SubSystem.h>

namespace dolfin
{

class BilinearForm;
class GenericMatrix;
class GenericVector;
class Mesh;
class SubDomain;

/// Common base class for boundary conditions

class BoundaryCondition : public MeshDependent
{

public:

  /// Apply boundary condition to linear system
  virtual void apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form) = 0;

  /// Apply boundary condition to linear system for a nonlinear problem
  virtual void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
                     BilinearForm const& form) = 0;

  /// Apply boundary condition to a subsystem of the linear system
  /// Implemented as changing the subsystem temporarily
  /// NOTE: the implementation of the boundary condition should not cache any
  ///       data structure assuming the subsystem is left unchanged
  virtual void apply(GenericMatrix& A, GenericVector& b,
                     BilinearForm const& form, SubSystem const sub_system);

  /// Apply boundary condition to a subsystem of the linear system for a nonlinear problem
  /// Implemented as changing the subsystem temporarily
  /// NOTE: the implementation of the boundary condition should not cache any
  ///       data structure assuming the subsystem is left unchanged
  virtual void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
                     BilinearForm const& form, SubSystem const sub_system);

  ///
  BoundaryCondition& operator()(Time const& t)
  {
    this->sync(t);
    return *this;
  }

  ///
  std::string const& type() const;

  ///
  Mesh& mesh() const;

  // Sub domain
  SubDomain const& sub_domain() const;

  ///
  SubSystem const& sub_system() const;

protected:

  /// Constructor based on a geometrical subdomain
  BoundaryCondition(std::string const& type, Mesh& mesh,
                    SubDomain const& sub_domain);

  /// Constructor on a geometrical subdomain for a given subspace
  BoundaryCondition(std::string const& type, Mesh& mesh,
                    SubDomain const& sub_domain, SubSystem const sub_system);

  ///
  virtual void sync(Time const& t) = 0;

public:

  /// Destructor
  virtual ~BoundaryCondition();

private:

  // Default constructor
  BoundaryCondition();

  // String identifier for the boundary condition type.
  std::string const type_;

  // Mesh
  Mesh& mesh_;

  // Sub domain
  SubDomain const& sub_domain_;

  // Sub system
  SubSystem sub_system_;

};

//-----------------------------------------------------------------------------
inline std::string const& BoundaryCondition::type() const
{
  return type_;
}
//-----------------------------------------------------------------------------
inline Mesh& BoundaryCondition::mesh() const
{
  return mesh_;
}
//-----------------------------------------------------------------------------
inline SubDomain const& BoundaryCondition::sub_domain() const
{
  return sub_domain_;
}
//-----------------------------------------------------------------------------
inline SubSystem const& BoundaryCondition::sub_system() const
{
  return sub_system_;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDARY_CONDITION_H */
