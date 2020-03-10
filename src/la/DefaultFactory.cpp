// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/DefaultFactory.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/la/JANPACKFactory.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
GenericMatrix* DefaultFactory::createMatrix() const
{
  return factory().createMatrix();
}
//-----------------------------------------------------------------------------
GenericVector* DefaultFactory::createVector() const
{
  return factory().createVector();
}
//-----------------------------------------------------------------------------
GenericSparsityPattern * DefaultFactory::createPattern() const
{
  return factory().createPattern();
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& DefaultFactory::factory()
{

  // Get backend from parameter system
  std::string backend = dolfin_get<std::string>("linear algebra backend");

#ifdef HAVE_PETSC
  if (backend == "PETSc")
  {
    return PETScFactory::instance();
  }
#endif

#ifdef HAVE_JANPACK
  if (backend == "JANPACK")
  {
    return JANPACKFactory::instance();
  }
#endif

  error("Linear algebra backend \"" + backend + "\" not available.");
  return factory(); // Never reached :P
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
