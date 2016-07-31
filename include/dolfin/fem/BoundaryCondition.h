// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007, 2008.
// Modified by Aurélien Larcher, 2014.
//
// This abstraction was rewritten to violate abstraction the least possible
// instead of "in every possible way".
// Additionally support for subsystems was added.
//
// First added:  2008-06-18
// Last changed: 2014-04-02

#ifndef __DOLFIN_BOUNDARY_CONDITION_H
#define __DOLFIN_BOUNDARY_CONDITION_H

#include <dolfin/mesh/MeshDependent.h>

#include <dolfin/fem/SubSystem.h>

#include <ufc.h>

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

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDARY_CONDITION_H */
