// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007, 2008.
//
// First added:  2008-06-18
// Last changed: 2007-12-09

#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "UFCMesh.h"
#include "DofMap.h"

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
class Form;

/// Common base class for boundary conditions

class BoundaryCondition
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
  virtual void apply(GenericMatrix& A, GenericVector& b, const Form& form) = 0;

  /// Apply boundary condition to linear system for a nonlinear problem
  virtual void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
                     const Form& form) = 0;

  ///
  std::string const& type() const;

  ///
  Mesh& mesh() const;

  ///
  bool const has_geometrical_sub_domain() const;

  ///
  SubDomain const& sub_domain() const;

  ///
  uint const& sub_domain_index() const;

  ///
  MeshFunction<uint> const& sub_domain_markers() const;

  ///
  SubSystem const& sub_system() const;

protected:

  // Mark sub domain with mesh function
  void init_markers(uint const& topological_dim);

  // Local data for application of boundary conditions
  class LocalData
  {

  public:

    // Constructor
    LocalData(Mesh& mesh, SubSystem const& sub_system, Form const& form);

    // Destructor
    ~LocalData();

    // UFC view of mesh
    UFCMesh ufc_mesh;

    // Finite element for sub system
    ufc::finite_element const * finite_element;

    // Dof map for sub system
    DofMap const * dof_map;

    // Offset for sub system
    uint offset;

    // Local data used to set boundary conditions
    real* w;
    uint* cell_dofs;
    uint* facet_dofs;
    real** coordinates;

  private:

      //
      bool const is_subspace_;

  };

  LocalData& updateLocalData(Form const& form) const;

private:

  // Default constructor
  BoundaryCondition();

  // String identifier for the boundary condition type.
  std::string const type_;

  // Mesh
  Mesh& mesh_;

  //
  mutable LocalData * data_;

  // Sub domain index
  uint const sub_domain_index_;

  // Is the subdomain geometrical
  bool const has_geometrical_sub_domain_;

  // Sub domain
  SubDomain const * geometrical_sub_domain_;

  // Sub domain markers
  MeshFunction<uint> * sub_domain_markers_;

  // True if sub domain markers are created locally
  bool const local_sub_domain_markers_;

  // Sub system
  SubSystem const sub_system_;

};

//--- INLINE ------------------------------------------------------------------

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
inline bool const BoundaryCondition::has_geometrical_sub_domain() const
{
  return has_geometrical_sub_domain_;
}

//-----------------------------------------------------------------------------
inline SubDomain const& BoundaryCondition::sub_domain() const
{
  return *geometrical_sub_domain_;
}

//-----------------------------------------------------------------------------
inline uint const& BoundaryCondition::sub_domain_index() const
{
  return sub_domain_index_;
}

//-----------------------------------------------------------------------------
inline MeshFunction<uint> const& BoundaryCondition::sub_domain_markers() const
{
  dolfin_assert(sub_domain_markers_);
  return *sub_domain_markers_;
}

//-----------------------------------------------------------------------------
inline SubSystem const& BoundaryCondition::sub_system() const
{
  return sub_system_;
}

//-----------------------------------------------------------------------------
inline BoundaryCondition::LocalData&
BoundaryCondition::updateLocalData(Form const& form) const
{
  if(data_ == NULL)
  {
    data_ = new LocalData(mesh_, sub_system_, form);
  }
  return *data_;
}

}

#endif
