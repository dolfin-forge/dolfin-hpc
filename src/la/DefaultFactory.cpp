// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-17
// Last changed: 2008-05-19

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/la/EpetraFactory.h>
#include <dolfin/la/JANPACKFactory.h>
#include <dolfin/la/DefaultFactory.h>

#ifndef NO_UBLAS
#include <dolfin/la/uBlasFactory.h>
#endif


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
#ifndef NO_UBLAS
  // Fallback
  std::string default_backend = "uBLAS";
  typedef uBlasFactory<> DefaultFactory;
#endif

  // Get backend from parameter system
  std::string backend = dolfin_get("linear algebra backend");

#ifndef NO_UBLAS
  // Choose backend
  if (backend == "uBLAS")
  {
    return uBlasFactory<>::instance();
  }
  else if (backend == "PETSc")
#else
  if (backend == "PETSc")
#endif
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

#ifndef NO_UBLAS 
  // Fallback
  message("Linear algebra backend \"" + backend + "\" not available, using " + default_backend + ".");
  return DefaultFactory::instance();
#endif
}
//-----------------------------------------------------------------------------
