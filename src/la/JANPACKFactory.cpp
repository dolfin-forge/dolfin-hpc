// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/JANPACKFactory.h>

#ifdef HAVE_JANPACK

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKMat* JANPACKFactory::createMatrix() const 
{ 
  JANPACKMat* jm = new JANPACKMat();
  return jm;
}
//-----------------------------------------------------------------------------
JANPACKVec* JANPACKFactory::createVector() const 
{   
  message("Create vector");
  JANPACKVec* jv = new JANPACKVec();
  return jv;
}
//-----------------------------------------------------------------------------

// Singleton instance
JANPACKFactory JANPACKFactory::factory;

#endif
