// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2009.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2007-05-01
// Last changed: 2010-05-09

#ifndef __SLIPBC_H
#define __SLIPBC_H

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

#include <ufc.h>

#include <set>

namespace dolfin
{

class BoundaryMesh;
class BoundaryNormal;
class DofMap;
class Form;
class Function;
class GenericMatrix;
class GenericVector;
class Mesh;
class SubDomain;

class SlipBC: public BoundaryCondition
{

public:

  /// Create boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain);

  /// Create boundary condition for sub domain given a boundary normal function
  SlipBC(BoundaryNormal& normal, SubDomain const& sub_domain);

  /// Create boundary condition for sub domain specified by index
  SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain);

  /// Create sub system boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain, SubSystem const& sub_system);

  /// Create sub system boundary condition for sub domain specified by index
  SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
         SubSystem const& sub_system);

  /// Destructor
  ~SlipBC();

  /// Access to boundary normals
  BoundaryNormal& normal();

  //--- INTERFACE -------------------------------------------------------------

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, Form const& form);

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, DofMap const& dof_map,
             ufc::form const& ufc_form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
             Form const& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
             DofMap const& dof_map, ufc::form const& ufc_form);

private:

  void applySlipBC(Matrix& A, Matrix& As, Vector&, Mesh& mesh, uint node,
                   Array<uint>& nodes);

  // Initialize sub domain markers
  void init(SubDomain const& sub_domain);

  void apply(GenericMatrix& A, GenericVector& b, DofMap const& dof_map,
             Form const& form);

  // Global mesh
  Mesh& mesh_;

  // Boundary mesh
  BoundaryMesh& boundary_;

  // Sub domain markers (if any)
  MeshFunction<uint>* sub_domains_;

  // The sub domain
  uint sub_domain_;

  // True if sub domain markers are created locally
  bool local_sub_domains_;

  // Sub system
  SubSystem sub_system_;

  // User defined sub domain
  SubDomain const * user_sub_domain_;

  // Node normal and tangents
  BoundaryNormal * const normal_;

  Matrix* As_;

  int N_local_;
  int N_offset_;
  std::set<uint> off_proc_rows_;

  // Local data structures for assembly
  uint const tdim_;
  Array<real> a[3];
  Array<uint> a_col_indices[3];
  uint a_ncols[3];
  real l[3];
  real basis_[3][3];
  uint maxcomp[3];

  static real const permutation_matrix_[3][3];

};

}

#endif
