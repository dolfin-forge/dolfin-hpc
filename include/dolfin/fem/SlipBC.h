// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2009.
// Modified by Aurélien Larcher, 2013-2014. (rewrite, extension to any element)
//
// First added:  2007-05-01
// Last changed: 2014-05-22

#ifndef __DOLFIN_SLIPBC_H
#define __DOLFIN_SLIPBC_H

#include "BoundaryCondition.h"
#include "NodeNormal.h"
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

#include <set>

namespace dolfin
{

class DofMap;
class Form;
class Function;
class Mesh;
class ScratchSpace;
class SubDomain;

class SlipBC : public BoundaryCondition
{
public:

  /// Create boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain);

  /// Create boundary condition for sub domain given normals
  SlipBC(Mesh& mesh, SubDomain const& sub_domain, NodeNormal& normals);

  /// Create sub system boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain, SubSystem const& sub_system);

  /// Destructor
  ~SlipBC();

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form);

  /// Apply boundary condition to non linear system
  void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
             BilinearForm const& form);

  BoundaryNormal& normal()
  {
    return *node_normal;
  }

private:

  inline void sync(Time const&) { /* No-op */ }

  void applySlipBC_P1(GenericMatrix& A, GenericVector& b,
                      BilinearForm const& form, ScratchSpace& scratch);

  void applySlipBC(GenericMatrix& A, GenericVector& b,
                   BilinearForm const& form, ScratchSpace& scratch);

  void applyNodeBC(GenericMatrix& A, GenericVector& b, Mesh const& mesh,
                   uint const node, Array<uint> const& udofs,
                   Array<uint> const& ndofs);

  // Initialize sub domain markers
  void init(SubDomain const& sub_domain);

  // The mesh
  Mesh& mesh;

  // Node normal and tangents
  NodeNormal * node_normal;
  bool node_normal_local;

  Matrix* As;

  // Local data structures for assembly
  std::set<uint> row_indices;
  Array<real> a[3]; // local lhs extracted from A
  Array<real> a_slip_row[3]; // local lhs row after slip enforcement
  Array<uint> a_col_indices[3]; // non-zero indices per row
  real l[3];  // local rhs extracted from b
  real l_slip[3];  // local rhs after slip enforcement
  real basis_[3][3]; // local basis (n, tau1, tau2)
  uint max[3]; // maximum component
  uint row[3]; // row reordering


};

} /* namespace dolfin */

#endif /* __DOLFIN_SLIPBC_H */
