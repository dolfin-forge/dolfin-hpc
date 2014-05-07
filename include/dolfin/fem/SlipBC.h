// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2010-05-09

#ifndef __SLIPBC_H
#define __SLIPBC_H

#include "BoundaryCondition.h"
#include "NodeNormal.h"
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

#include <set>

namespace dolfin
{
class DofMap;
class Function;
class Mesh;
class SubDomain;
class Form;

class SlipBC : public BoundaryCondition
{
public:

  /// Create boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain);

  /// Create boundary condition for sub domain given normals
  SlipBC(Mesh& mesh, SubDomain const& sub_domain, NodeNormal& normals);

  /// Create boundary condition for sub domain specified by index
  SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain);

  /// Create sub system boundary condition for sub domain
  SlipBC(Mesh& mesh, SubDomain const& sub_domain, const SubSystem& sub_system);

  /// Create sub system boundary condition for sub domain specified by index
  SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain,
         const SubSystem& sub_system);

  /// Destructor
  ~SlipBC();

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form);

  /// Apply boundary condition to non linear system
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             const BilinearForm& form);

  BoundaryNormal& normal()
  {
    return *node_normal;
  }

private:

  void applySlipBC(Matrix& A, Matrix& As, Vector&, Mesh& mesh, uint node,
                   Array<uint>& nodes);

  // Do: A(row,col) = value   using setblock not setvalue
  inline void Aset(Matrix& A, uint row, uint col, real value)
  {
    A.set(&value, 1, &row, 1, &col);
  }

  // Do: b(row) = value   using setblock not setvalue
  inline void bset(Vector& b, uint row, real value)
  {
    b.set(&value, 1, &row);
  }

  // Initialize sub domain markers
  void init(SubDomain const& sub_domain);

  // The mesh
  Mesh& mesh;

  // Node normal and tangents
  NodeNormal * node_normal;
  bool node_normal_local;

  Matrix* As;
  std::set<uint> off_proc_rows;

  real *row_block;
  real *zero_block;
  uint *a1_indices_array;
};

}

#endif
