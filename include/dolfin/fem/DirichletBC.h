// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Kristian Oelgaard, 2007
//
// First added:  2007-04-10
// Last changed: 2008-05-22
//

#ifndef __DOLFIN_DIRICHLET_BC_H
#define __DOLFIN_DIRICHLET_BC_H

#include <dolfin/fem/BoundaryCondition.h>

#include <dolfin/common/types.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>

namespace dolfin
{

class DofMap;
class Coefficient;
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
 *  For mixed systems (vector-valued and mixed elements), an
 *  optional set of parameters may be used to specify for which sub
 *  system the boundary condition should be specified.
 */

class DirichletBC: public BoundaryCondition
{
public:

  /// Create boundary condition for sub domain
  DirichletBC(Coefficient& g, Mesh& mesh, const SubDomain& sub_domain,
              BCMethod method = topological);

  /// Create sub system boundary condition for sub domain
  DirichletBC(Coefficient& g, Mesh& mesh, const SubDomain& sub_domain,
              const SubSystem& sub_system, BCMethod method = topological);

  /// Destructor
  ~DirichletBC();

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
             BilinearForm const& form);

private:

  /// Apply boundary conditions
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector* x,
             BilinearForm const& form);

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
  Coefficient& g_;

  // Search method
  BCMethod method_;

  // Boundary facets, stored as pairs (cell, local facet number)
  Array<uint> * entities_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_DIRICHLET_BC_H */
