// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007
//
// First added:  2007-07-08
// Last changed: 2007-12-08

#ifndef __DOLFIN_PERIODIC_BC_H
#define __DOLFIN_PERIODIC_BC_H

#include <dolfin/common/types.h>
#include "SubSystem.h"
#include "BoundaryCondition.h"

namespace dolfin
{

class DofMap;
class Mesh;
class SubDomain;
class Form;
class GenericMatrix;
class GenericVector;
class PeriodicSubDomain;

/// This class specifies the interface for setting periodic boundary
/// conditions for partial differential equations,
///
///    u(x) = u(F^{-1}(x)) on G,
///    u(x) = u(F(x))      on H,
///
/// where F : H --> G is a map from a subdomain H to a subdomain G.
///
/// A PeriodicBC is specified by a Mesh and a SubDomain. The given
/// subdomain must overload both the inside() function, which
/// specifies the points of G, and the map() function, which
/// specifies the map from the points of H to the points of G.
///
/// For mixed systems (vector-valued and mixed elements), an
/// optional set of parameters may be used to specify for which sub
/// system the boundary condition should be specified.

class PeriodicBC: public BoundaryCondition
{
public:

  /// Create periodic boundary condition for sub domain
  PeriodicBC(Mesh& mesh, PeriodicSubDomain const& sub_domain);

  /// Create sub system boundary condition for sub domain
  PeriodicBC(Mesh& mesh, PeriodicSubDomain const& sub_domain,
             SubSystem const& sub_system);

  /// Destructor
  ~PeriodicBC();

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form);

  /// Apply boundary condition to linear system for a nonlinear problem (not implemented)
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             BilinearForm const& form);

private:

};

} /* namespace dolfin */

#endif /* __DOLFIN_PERIODIC_BC_H */
