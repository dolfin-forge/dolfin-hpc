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

#include <dolfin/common/types.h>
#include <dolfin/fem/SubSystem.h>

#include <ufc.h>

namespace dolfin
{

class DofMap;
class GenericMatrix;
class GenericVector;
class SubDomain;
class Mesh;
class BilinearForm;
template<class T> class MeshFunction;

/// Common base class for boundary conditions

class BoundaryCondition : public MeshDependent
{
public:

  /// Constructor based on a geometrical subdomain
  BoundaryCondition(std::string const& type, Mesh& mesh,
                    SubDomain const& sub_domain);

  /// Constructor based on boundary sub domain markers
  BoundaryCondition(std::string const& type, MeshFunction<uint>& sub_domains,
                    uint sub_domain);

  /// Constructor on a geometrical subdomain for a given subspace
  BoundaryCondition(std::string const& type, Mesh& mesh,
                    SubDomain const& sub_domain, SubSystem const sub_system);

  /// Constructor on boundary sub domain markers for a given subspace
  BoundaryCondition(std::string const& type, MeshFunction<uint>& sub_domains,
                    uint sub_domain, SubSystem const sub_system);

  /// Destructor
  virtual ~BoundaryCondition();

  /// Apply boundary condition to linear system
  virtual void apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form) = 0;

  /// Apply boundary condition to linear system for a nonlinear problem
  virtual void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                     const BilinearForm& form) = 0;

  /// Apply boundary condition to a subsystem of the linear system
  /// Implemented as changing the subsystem temporarily
  /// NOTE: the implementation of the boundary condition should not cache any
  ///       data structure assuming the subsystem is left unchanged
  virtual void apply(GenericMatrix& A, GenericVector& b,
                     const BilinearForm& form, SubSystem const sub_system);

  /// Apply boundary condition to a subsystem of the linear system for a nonlinear problem
  /// Implemented as changing the subsystem temporarily
  /// NOTE: the implementation of the boundary condition should not cache any
  ///       data structure assuming the subsystem is left unchanged
  virtual void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                     const BilinearForm& form, SubSystem const sub_system);

  ///
  std::string const& type() const;

  ///
  Mesh& mesh() const;

  ///
  SubSystem const& sub_system() const;

protected:

  ///
  bool has_geometrical_sub_domain() const;

  ///
  SubDomain const& sub_domain() const;

  ///
  uint const& sub_domain_index() const;

  ///
  MeshFunction<uint> const& sub_domain_markers() const;

  /// Mark sub domain with mesh function defined for given mesh entity type
  /// Convention: 0 for matching entities, 1 elsewhere
  void init_markers(uint const entity_dim);

private:

  // Default constructor
  BoundaryCondition();

  // String identifier for the boundary condition type.
  std::string const type_;

  // Mesh
  Mesh& mesh_;

  // Sub domain index
  uint const sub_domain_index_;

  // Is the subdomain geometrical
  bool const has_geometric_sub_domain_;

  // Sub domain
  SubDomain const * geometric_sub_domain_;

  // Sub domain markers
  MeshFunction<uint> * sub_domain_markers_;

  // True if sub domain markers are created locally
  bool const local_sub_domain_markers_;

  // Sub system
  SubSystem sub_system_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDARY_CONDITION_H */
