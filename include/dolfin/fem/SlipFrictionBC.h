// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-07-15 (merged from branch larcher)
// Last changed: 2013-07-15

#ifndef __DOLFIN_SLIP_FRICTION_BC_H
#define __DOLFIN_SLIP_FRICTION_BC_H

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/SlipBC.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>

#include <ufc.h>

namespace dolfin
{

class Form;
class Function;
class GenericMatrix;
class SubDomain;

class SlipFrictionBC: public BoundaryCondition
{
public:
  /// Create boundary condition for sub domain
  SlipFrictionBC(Mesh& mesh, SubDomain const& sub_domain, real beta);

  /// Create sub system boundary condition for sub domain
  SlipFrictionBC(Mesh& mesh, SubDomain const& sub_domain,
                 SubSystem const& sub_system, real beta);

  /// Create boundary condition for sub domain specified by index
//  SlipFrictionBC(MeshFunction<uint>& sub_domains, uint sub_domain, real beta);

  /// Create sub system boundary condition for sub domain specified by index
//  SlipFrictionBC(MeshFunction<uint>& sub_domains, uint sub_domain,
//                 SubSystem const& sub_system, real beta);

/// Destructor
  ~SlipFrictionBC();

  BoundaryNormal& normal();

  //--- INTERFACE -------------------------------------------------------------

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, const BilinearForm& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             const BilinearForm& form);

  real beta() const;

  Function& friction();

private:

  SlipBC slipbc_;

  real beta_;
  IndicatorExpression expr_;
  Function Fbeta_;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline void SlipFrictionBC::apply(GenericMatrix& A, GenericVector& b,
                                  BilinearForm const& form)
{
  slipbc_.apply(A, b, form);
}

//-----------------------------------------------------------------------------
inline void SlipFrictionBC::apply(GenericMatrix& A, GenericVector& b,
                                  const GenericVector& x, BilinearForm const& form)
{
  slipbc_.apply(A, b, x, form);
}

//-----------------------------------------------------------------------------
inline real SlipFrictionBC::beta() const
{
  return beta_;
}

//-----------------------------------------------------------------------------
inline Function& SlipFrictionBC::friction()
{
  return Fbeta_;
}

} /* namespace dolfin */

#endif /* __DOLFIN_SLIP_FRICTION_BC_H */

