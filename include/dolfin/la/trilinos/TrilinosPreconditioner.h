// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_PRECONDITIONER_H
#define __DOLFIN_TRILINOS_PRECONDITIONER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosObject.h>

namespace dolfin
{

namespace trilinos
{

class Vector;
enum class PreconditionerType;

/// This class specifies the interface for user-defined Krylov method
/// Preconditioners. A user wishing to implement their own Preconditioner needs
/// only supply a function that approximately solves the linear system given a
/// right-hand side.

class Preconditioner : public Object
{
public:
  /// Constructor
  Preconditioner() = default;

  /// Destructor
  virtual ~Preconditioner() = default;

  // static void setup( const KSP ksp, Preconditioner & pc );

  /// Solve linear system approximately for given right-hand side b
  virtual void solve( trilinos::Vector & x, trilinos::Vector const & b ) = 0;

  /// Friends
  friend class KrylovSolver;

protected:
  // PC petscpc;

private:
  // static auto PCApply( PC pc, Vec x, Vec y ) -> int;
  // static auto PCCreate( PC pc ) -> int;

  /// Return trilinos Preconditioner type
  // static auto getType( PreconditionerType pc ) -> PCType;
};

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_PRECONDITIONER_H
