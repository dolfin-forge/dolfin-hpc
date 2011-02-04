// Copyright (C) 2005-2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-12-06
// Last changed: 2007-12-07

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/PETScFactory.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
SparsityPattern* PETScFactory::createPattern() const 
{
  return new SparsityPattern(); 
}
//-----------------------------------------------------------------------------
PETScMatrix* PETScFactory::createMatrix() const 
{ 
  PETScMatrix* pm = new PETScMatrix();
  return pm;
}
//-----------------------------------------------------------------------------
PETScVector* PETScFactory:: createVector() const 
{ 
  return new PETScVector(); 
}
//-----------------------------------------------------------------------------

// Singleton instance
PETScFactory PETScFactory::factory;

#endif
