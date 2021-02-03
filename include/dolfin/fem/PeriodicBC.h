// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PERIODIC_BC_H
#define __DOLFIN_PERIODIC_BC_H

#include <dolfin/common/types.h>
#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/SubSystem.h>

namespace dolfin
{

class Mesh;
class BilinearForm;
class GenericMatrix;
class GenericVector;
class PeriodicSubDomain;

//-----------------------------------------------------------------------------

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

class PeriodicBC : public BoundaryCondition
{
public:
  /// Create periodic boundary condition for sub domain
  PeriodicBC( Mesh & mesh, PeriodicSubDomain const & sub_domain );

  /// Create sub system boundary condition for sub domain
  PeriodicBC( Mesh &                    mesh,
              PeriodicSubDomain const & sub_domain,
              SubSystem const &         sub_system );

  /// Destructor
  ~PeriodicBC() override = default;

  /// Apply boundary condition to linear system
  void apply( GenericMatrix &      A,
              GenericVector &      b,
              BilinearForm const & form ) override;

  /// Apply boundary condition to linear system for a nonlinear problem (not
  /// implemented)
  void apply( GenericMatrix &       A,
              GenericVector &       b,
              const GenericVector & x,
              BilinearForm const &  form ) override;
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_PERIODIC_BC_H */
