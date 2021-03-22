// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/DefaultFactory.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/janpack/JANPACKFactory.h>
#include <dolfin/la/petsc/PETScFactory.h>
#include <dolfin/la/trilinos/TrilinosFactory.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

auto DefaultFactory::createMatrix() const -> GenericMatrix *
{
  return factory().createMatrix();
}

//-----------------------------------------------------------------------------

auto DefaultFactory::createVector() const -> GenericVector *
{
  return factory().createVector();
}

//-----------------------------------------------------------------------------

auto DefaultFactory::createPattern() const -> GenericSparsityPattern *
{
  return factory().createPattern();
}

//-----------------------------------------------------------------------------

auto DefaultFactory::factory() -> LinearAlgebraFactory &
{

  // Get backend from parameter system
  std::string backend = dolfin_get< std::string >( "linear algebra backend" );

#ifdef HAVE_PETSC
  if ( backend == "PETSc" )
  {
    return PETScFactory::instance();
  }
#endif

#ifdef HAVE_TRILINOS
  if ( backend == "Trilinos" )
  {
    return trilinos::Factory::instance();
  }
#endif

#ifdef HAVE_JANPACK
  if ( backend == "JANPACK" )
  {
    return JANPACKFactory::instance();
  }
#endif

  error( "Linear algebra backend \"" + backend + "\" not available." );
  return factory(); // Never reached :P
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
