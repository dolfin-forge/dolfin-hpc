// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2010-2011.
//
// First added:  2008-05-17
// Last changed: 2011-06-10

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/la/EpetraFactory.h>
#include <dolfin/la/JANPACKFactory.h>
#include <dolfin/la/DefaultFactory.h>


using namespace dolfin;

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
LinearAlgebraFactory& DefaultFactory::factory() const
{

  // Get backend from parameter system
  std::string backend = dolfin_get("linear algebra backend");

  if (backend == "PETSc")
  {
#ifdef HAVE_PETSC
    return PETScFactory::instance();
#endif
  }
  else if (backend == "Epetra")
  {
#ifdef HAVE_TRILINOS
    return EpetraFactory::instance();
#endif
  }
  else if (backend == "JANPACK")
  {
#ifdef HAVE_JANPACK
    return JANPACKFactory::instance();
#endif
  }

  error("Linear algebra backend \"" + backend + "\" not available.");

}
//-----------------------------------------------------------------------------
