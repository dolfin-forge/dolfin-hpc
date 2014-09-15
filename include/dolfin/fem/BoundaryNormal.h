// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#ifndef __BOUNDARYNORMAL_H
#define __BOUNDARYNORMAL_H

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/function/Function.h>

namespace dolfin
{

class BoundaryMesh;
class FiniteElementSpace;
class Mesh;

/**
 *  DOCUMENTATION:
 *
 *  @class  BoundaryNormal
 *
 *  @brief  Provides an abstraction for a normal function defined at the
 *          boundary of the computational domain, and the tangential vector(s).
 */

class BoundaryNormal
{
public:

  /// Return boundary mesh
  BoundaryMesh& boundary();

  /// Destructor
  virtual ~BoundaryNormal();

  /// Return global mesh
  Mesh& mesh();

  /// Return the basis (n, tau, 0) in 2d or (n, tau1, tau2) in 3d
  Array<Function>& basis();

  /// Return the node type
  Function& node_type();

  /// Initialization of basis functions and node type for given space
  void init(FiniteElementSpace const& space);

  /// Write orthonormal basis and node type to file
  void write(std::string const& filename);

  ///
  virtual void compute() = 0;

protected:

  ///
  BoundaryNormal(Mesh& mesh);

  ///
  explicit BoundaryNormal(BoundaryMesh& boundary);

private:

  Mesh& mesh_;
  BoundaryMesh * const boundary_;
  bool const local_boundary_;
  mutable Array<Function> basis_;
  mutable Function node_type_;

};

}
#endif

