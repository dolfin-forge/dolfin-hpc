// Copyright (C) 2005-2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/PETScFactory.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
auto PETScFactory::createPattern() const -> SparsityPattern* 
{
  return new SparsityPattern(); 
}
//-----------------------------------------------------------------------------
auto PETScFactory::createMatrix() const -> PETScMatrix* 
{ 
  PETScMatrix* pm = new PETScMatrix();
  return pm;
}
//-----------------------------------------------------------------------------
auto PETScFactory:: createVector() const -> PETScVector* 
{ 
  return new PETScVector(); 
}
//-----------------------------------------------------------------------------

// Singleton instance
PETScFactory PETScFactory::factory;

#endif
