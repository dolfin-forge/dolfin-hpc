// Copyright (C) 2005-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-09-02
// Last changed: 2008-06-23

#include <dolfin/evolution/TimeDependent.h>

#include <dolfin/evolution/Time.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
TimeDependent::TimeDependent(Time const& time) :
    t_(time.clock())
{
  // Do nothing
}
//-----------------------------------------------------------------------------
TimeDependent::TimeDependent(real const& t) :
    t_(t)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
TimeDependent::~TimeDependent()
{
  // Do nothing
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
