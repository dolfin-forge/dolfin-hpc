// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-07-15 (merged from branch larcher)
// Last changed: 2013-07-15

#ifndef __DOLFIN_SLIP_FRICTION_BC_H
#define __DOLFIN_SLIP_FRICTION_BC_H

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/SlipBC.h>

#include <ufc.h>

namespace dolfin
{

class Form;
class Function;
class GenericMatrix;
class SubDomain;

class SlipFrictionBC : public BoundaryCondition
{
public:

  /// Create boundary condition for sub domain
  SlipFrictionBC(Coefficient& beta, Mesh& mesh, SubDomain const& sub_domain);

  /// Create boundary condition for sub domain specified by index
  SlipFrictionBC(Coefficient& beta, MeshFunction<uint>& sub_domains,
                 uint sub_domain);

  /// Create sub system boundary condition for sub domain
  SlipFrictionBC(Coefficient& beta, Mesh& mesh, SubDomain const& sub_domain,
                 SubSystem const& sub_system);

  /// Create sub system boundary condition for sub domain specified by index
  SlipFrictionBC(Coefficient& beta, MeshFunction<uint>& sub_domains,
                 uint sub_domain, SubSystem const& sub_system);

  /// Destructor
  ~SlipFrictionBC();

  ///
  BoundaryNormal& normal();

  ///
  Coefficient& friction() const;

  //--- INTERFACE -------------------------------------------------------------

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, BilinearForm const& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, GenericVector const& x,
             BilinearForm const& form);

private:

  SlipBC slipbc_;
  Coefficient& beta_;

};

//--- INLINES -----------------------------------------------------------------

inline void SlipFrictionBC::apply(GenericMatrix& A, GenericVector& b,
                                  BilinearForm const& form)
{
  slipbc_.apply(A, b, form);
}

//-----------------------------------------------------------------------------
inline void SlipFrictionBC::apply(GenericMatrix& A, GenericVector& b,
                                  GenericVector const& x,
                                  BilinearForm const& form)
{
  slipbc_.apply(A, b, x, form);
}

//-----------------------------------------------------------------------------
inline Coefficient& SlipFrictionBC::friction() const
{
  return beta_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_SLIP_FRICTION_BC_H */

