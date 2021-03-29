// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosPreconditioner.h>

#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/trilinos/TrilinosVector.h>
#include <dolfin/la/trilinos/TrilinosMatrix.h>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P, std::string pc_str ) -> void
{
  init( P, pc_type( pc_str ) );
}

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P, PreconditionerType pc_type ) -> void
{
  if ( pc_type == PreconditionerType::cheb )
  {
    _name = "CHEBYSHEV";
    // FIXME set parameters?
  }
  else if ( pc_type == PreconditionerType::ilut )
  {
    _name = "ILUT";
    // FIXME set parameters?
  }
  else if ( pc_type == PreconditionerType::riluk )
  {
    _name = "RILUK";
    // FIXME set parameters?
  }
  else if ( pc_type == PreconditionerType::relax )
  {
    _name = "RELAXATION";
    // FIXME set parameters?
  }
  else
  {
    warning( "Undefined preconditioner          "
             "Fallback to default preconditioner" );

    _name = "RELAXATION";
    // FIXME set parameters?
  }

  init( P );
}

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P ) -> void
{
  Ifpack2::Factory pc;

  pc_ = pc.create( _name, P.mat().getConst() );
  Teuchos::RCP< Teuchos::ParameterList > plist = Teuchos::parameterList();
  // FIXME set some parameters here?
  pc_->setParameters( *plist );
  pc_->initialize();
  pc_->compute();
}

//-----------------------------------------------------------------------------

auto Preconditioner::set( KrylovSolver & solver ) -> void
{
  // solver._problem->setRightPrec( pc_ );
}

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
