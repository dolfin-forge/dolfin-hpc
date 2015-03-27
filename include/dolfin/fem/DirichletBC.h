// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Kristian Oelgaard, 2007
//
// First added:  2007-04-10
// Last changed: 2008-05-22
//
// FIXME: This class needs some cleanup, in particular collecting
// FIXME: all data from different representations into a common
// FIXME: data structure (perhaps an Array<uint> with facet indices).

#ifndef __DIRICHLET_BC_H
#define __DIRICHLET_BC_H

#include <dolfin/common/types.h>
#include "SubSystem.h"
#include <dolfin/mesh/MeshFunction.h>
#include "BoundaryCondition.h"
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>

namespace dolfin
{

class DofMap;
class Function;
class Mesh;
class FiniteElementSpace;
class Form;
class GenericMatrix;
class GenericVector;

/// The BCMethod variable may be used to specify the type of method
/// used to identify degrees of freedom on the boundary. Available
/// methods are: topological approach (default), geometric approach,
/// and pointwise approach. The topological approach is faster,
/// but will only identify degrees of freedom that are located on a
/// facet that is entirely on the boundary. In particular, the
/// topological approach will not identify degrees of freedom
/// for discontinuous elements (which are all internal to the cell).
/// A remedy for this is to use the geometric approach. To apply
/// pointwise boundary conditions e.g. pointloads, one will have to
/// use the pointwise approach which in turn is the slowest of the
/// three possible methods.

enum BCMethod
{
  topological, geometric, pointwise
};

/**
 *  @class  DirichletBc
 *
 *  @brief  This class specifies the interface for setting (strong)
 *  Dirichlet boundary conditions for partial differential equations,
 *
 *     u = g on G,
 *
 *  where u is the solution to be computed, g is a function
 *  and G is a sub domain of the mesh.
 *
 *  A DirichletBC is specified by the Function g, the Mesh,
 *  and boundary indicators on (a subset of) the mesh boundary.
 *
 *  The boundary indicators may be specified in a number of
 *  different ways.
 *
 *  The simplest approach is to specify the boundary by a SubDomain
 *  object, using the inside() function to specify on which facets
 *  the boundary conditions should be applied.
 *
 *  Alternatively, the boundary may be specified by a MeshFunction
 *  labeling all mesh facets together with a number that specifies
 *  which facets should be included in the boundary.
 *
 *  The third option is to attach the boundary information to the
 *  mesh. This is handled automatically when exporting a mesh from
 *  for example VMTK.
 *
 *  For mixed systems (vector-valued and mixed elements), an
 *  optional set of parameters may be used to specify for which sub
 *  system the boundary condition should be specified.
 */

class DirichletBC: public BoundaryCondition
{
public:

  /// Create boundary condition for sub domain
  DirichletBC(Function& g, Mesh& mesh, const SubDomain& sub_domain,
              BCMethod method = topological);

  /// Create boundary condition for sub domain specified by index
  DirichletBC(Function& g, MeshFunction<uint>& sub_domains, uint sub_domain,
              BCMethod method = topological);

  /// Create sub system boundary condition for sub domain
  DirichletBC(Function& g, Mesh& mesh, const SubDomain& sub_domain,
              const SubSystem& sub_system, BCMethod method = topological);

  /// Create sub system boundary condition for sub domain specified by index
  DirichletBC(Function& g, MeshFunction<uint>& sub_domains, uint sub_domain,
              const SubSystem& sub_system, BCMethod method = topological);

  /// Destructor
  ~DirichletBC();

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             const BilinearForm& form);

private:

  /// Apply boundary conditions
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector* x,
             const BilinearForm& form);

  // Initialize sub domain markers from sub domain
  void initFromSubDomain(SubDomain const& sub_domain);

  // Initialize sub domain markers from MeshFunction
  void initFromMeshFunction(MeshFunction<uint> const& sub_domains,
                            uint sub_domain);

  // Compute boundary values for facet (topological approach)
  void computeBCTopological(_map<uint, real>& boundary_values,
                            FiniteElementSpace const& space,
                            SubSystem const& sub_system);

  // Compute boundary values for facet (geometrical approach)
  void computeBCGeometric(_map<uint, real>& boundary_values,
                          FiniteElementSpace const& space,
                          SubSystem const& sub_system);

  // Compute boundary values for facet (pointwise approach)
  void computeBCPointwise(_map<uint, real>& boundary_values,
                          FiniteElementSpace const& space,
                          SubSystem const& sub_system);

  // The function
  Function& g_;

  // Search method
  BCMethod method_;

  // Boundary facets, stored as pairs (cell, local facet number)
  std::vector< std::pair<uint, uint> > facets_;

};

}

#endif
