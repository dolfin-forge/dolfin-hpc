// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-07-04
// Last changed: 2006-07-04

#include <dolfin/config/dolfin_config.h>

#ifndef NO_UBLAS
#include <dolfin/la/uBlasVector.h>
#include <dolfin/la/uBlasDummyPreconditioner.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
uBlasDummyPreconditioner::uBlasDummyPreconditioner() : uBlasPreconditioner()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
uBlasDummyPreconditioner::~uBlasDummyPreconditioner()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void uBlasDummyPreconditioner::solve(uBlasVector& x, const uBlasVector& b) const
{
  x.vec().assign(b.vec());
}
//-----------------------------------------------------------------------------
#endif
