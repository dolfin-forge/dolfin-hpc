// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_PRECONDITIONER_H
#define __DOLFIN_TRILINOS_PRECONDITIONER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosObject.h>
#include <dolfin/la/trilinos/TrilinosVector.h>

#include <Ifpack2_Factory.hpp>

#include <string>

namespace dolfin
{

enum class PreconditionerType;

namespace trilinos
{

class Matrix;
class KrylovSolver;

// Class implementing an interface to trilinos (Ifpack2) preconditioners

class Preconditioner : public Object
{
public:
  /// Preconditioner Type
  using PCType = Ifpack2::Preconditioner< real, int, size_t, Vector::TPNode >;

public:
  /// Constructor
  Preconditioner() = default;

  /// Destructor
  virtual ~Preconditioner() = default;

  /// Set the preconditioner type on a solver
  void set( trilinos::KrylovSolver & solver );

  /// Initialise preconditioner based on Operator P and preconditioner string
  auto init( Matrix const & P, std::string pc_str ) -> void;

  /// Initialise preconditioner based on Operator P and preconditioner type
  auto init( Matrix const & P, PreconditionerType pc_type ) -> void;

  /// Initialise preconditioner based on Operator P
  auto init( Matrix const & P ) -> void;

private:
  // name of the preconditioner
  std::string _name;

  // Ifpack2 preconditioner, to be constructed from a Tpetra Operator or Matrix
  Teuchos::RCP< PCType > pc_;
};

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_PRECONDITIONER_H
