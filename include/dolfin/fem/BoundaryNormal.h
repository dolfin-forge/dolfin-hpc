// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BOUNDARYNORMAL_H
#define __DOLFIN_BOUNDARYNORMAL_H

#include <dolfin/common/types.h>
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
  auto boundary() -> BoundaryMesh &;

  /// Destructor
  virtual ~BoundaryNormal();

  /// Return global mesh
  auto mesh() -> Mesh &;

  /// Return the basis (n, tau, 0) in 2d or (n, tau1, tau2) in 3d
  auto basis() -> std::vector< Function > &;

  /// Initialization of basis functions and node type for given space
  auto init( FiniteElementSpace const & space ) -> void;

  /// Write orthonormal basis and node type to file
  auto write( std::string const & filename ) -> void;

  ///
  virtual auto compute() -> void = 0;

protected:
  ///
  BoundaryNormal( Mesh & mesh );

  ///
  explicit BoundaryNormal( BoundaryMesh & boundary );

private:
  Mesh &                  mesh_;
  BoundaryMesh * const    boundary_;
  bool const              local_boundary_;
  std::vector< Function > basis_;
  Function                node_type_;
};

//-----------------------------------------------------------------------------
inline auto BoundaryNormal::mesh() -> Mesh &
{
  return mesh_;
}
//-----------------------------------------------------------------------------
inline auto BoundaryNormal::boundary() -> BoundaryMesh &
{
  return *boundary_;
}
//-----------------------------------------------------------------------------
inline auto BoundaryNormal::basis() -> std::vector< Function > &
{
  return basis_;
}

}
#endif
